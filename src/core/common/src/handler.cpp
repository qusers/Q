#include "handler.h"

#include <chrono>

#include "native_out.h"
#include "offdiag.h"
#include "std_output.h"

void Handler::run_iteration(int iteration) {
    reset_energies();

    // 3. nonbonded forces
    calc_nonbonded_forces();

    // 2. bonded forces and some constraints calculations
    calc_internal_forces(iteration);

    prepare_offdiags(iteration);

    // 5. leapfrog integration
    calc_leapfrog();

    // 6. update temperature after integration.
    // todo: Why calculate temperature again here?
    calc_temperature();

    // 8. print output to files
    print_outputs(iteration);
}

void Handler::initialize(const CommandInfo& command) {
    ctx.command_info = command;
    const ParseResult parsed = ctx.get_parse_result();
    initialize_system(parsed);
    configure(stage_from(parsed));  // fresh start: velocities are generated here
}

void Handler::initialize_system(const ParseResult& parsed) {
    if (system_loaded_) {
        throw InputError("Handler::initialize_system called twice; create a new Handler instead.");
    }

    system_ = parsed;
    ctx.init(system_);  // topology, FEP, masses, shell, LJ matrix, start coords/velocities

    // Depends only on System fields (cut-offs, LRF, charge groups, pair masks): build once.
    nonbonded_force_ = create_nonbonded_force_backend();
    nonbonded_force_->init(ctx, system_);

    system_loaded_ = true;
    configured_ = false;
    // Restart: coords/velocities come from the restart file. As in Fortran Q, a given
    // initial_temperature (with seed > 0) still draws new velocities on the first configure().
    const bool new_velocities = ctx.md.generate_velocities && ctx.md.random_seed > 0;
    state_ready_ = !ctx.fresh_start && !new_velocities;
}

void Handler::rebuild_constraints_and_bonded() {
    constraint_force_ = create_constraint_force_backend();
    constraint_force_->init(ctx, system_);

    // Bonded skips constrained bonds/angles, so it must follow the constraint set.
    bonded_force_ = create_bonded_force_backend();
    bonded_force_->init(ctx, system_, constraint_force_->data(), nonbonded_force_->special_pairs());
}

void Handler::configure(const StageSettings& stage) {
    if (!system_loaded_) {
        throw InputError("Handler::configure called before initialize_system.");
    }

    // Read theta_corr from the CURRENT water boundary before it is replaced.
    std::vector<double> theta = configured_ ? current_theta_corr() : system_.restart_theta_corr;

    ctx.apply_stage(stage);  // md stage fields, dt, lambdas upload, native_output

    const std::array<bool, 3> flags{
        ctx.md.constraint_solvent,
        ctx.md.constraint_solute,
        ctx.md.constraint_hydrogens,
    };
    if (!configured_ || flags != constraint_flags_) {
        rebuild_constraints_and_bonded();
        constraint_flags_ = flags;
    }

    restraint_force_ = create_restraint_force_backend();
    restraint_force_->init(ctx, stage.restrspos, stage.restrseqs, stage.restrdists,
                           stage.restrangs, stage.restrwalls);

    // Order matters: water boundary and integrator keep raw pointers to these.
    temperature_ = create_temperature_backend();
    temperature_->init(ctx, *constraint_force_);

    water_boundary_force_ = create_water_boundary_force_backend();
    water_boundary_force_->init(ctx, *temperature_, system_, theta);  // cstb depends on lambdas

    integrator_ = create_integrator_backend();
    integrator_->init(ctx, *constraint_force_, *temperature_);

    configured_ = true;

    // First stage of a fresh start: generate velocities once. Needs the constraints and masses built above.
    if (!state_ready_) {
        start_fresh(ctx.md.random_seed, ctx.md.initial_temperature);
    }
}

void Handler::start_fresh(int seed, double initial_temperature) {
    if (!configured_) {
        throw InputError("Handler::start_fresh called before configure.");
    }
    if (seed <= 0) {
        throw InputError("start_fresh requires a positive random seed.");
    }

    ctx.generate_velocities(seed, initial_temperature);  // reuses the existing buffer and uploads
    constraint_force_->initial_constraint(ctx);
    stop_cm_translation();
    state_ready_ = true;
}

void Handler::set_state(const SystemState& s) {
    if (!system_loaded_) {
        throw InputError("Handler::set_state called before load_system.");
    }
    ctx.set_state(s);  // validates sizes, fills cpu buffers, uploads coords + velocities

    // theta_corr lives in the water boundary backend: re-apply it there.
    if (configured_) {
        water_boundary_force_ = create_water_boundary_force_backend();
        water_boundary_force_->init(ctx, *temperature_, system_, s.theta_corr);
        integrator_ = create_integrator_backend();  // not strictly needed, keeps pointers fresh
        integrator_->init(ctx, *constraint_force_, *temperature_);
    } else {
        system_.restart_theta_corr = s.theta_corr;  // picked up by the first configure()
    }
    state_ready_ = true;
    broken_ = false;
}

std::vector<double> Handler::current_theta_corr() {
    std::vector<double> theta;
    if (!water_boundary_force_) return theta;

    water_boundary_force_->sync_for_output(ctx);  // copies backend wshells into ctx.wshells
    theta.reserve(ctx.wshells.size());
    for (const auto& shell : ctx.wshells) {
        theta.push_back(shell.theta_corr);
    }
    return theta;
}

SystemState Handler::get_state() {
    if (!configured_) {
        throw InputError("Handler::get_state called before configure.");
    }
    if (ctx.command_info.requested_gpu) {
        ctx.coords->download();
        ctx.velocities->download();
    }
    SystemState s;
    s.coords.assign(ctx.coords->cpu_data_p, ctx.coords->cpu_data_p + ctx.n_atoms);
    s.velocities.assign(ctx.velocities->cpu_data_p, ctx.velocities->cpu_data_p + ctx.n_atoms);
    s.theta_corr = current_theta_corr();
    return s;
}

void Handler::stop_cm_translation() {
    auto& velocities = ctx.velocities->cpu_data_p;
    auto& masses = ctx.masses->cpu_data_p;
    double total_mass = 0;
    coord_t vcm = {};

    for (int ai = 0; ai < ctx.n_atoms; ai++) {
        const double rmass = masses[ai];
        total_mass += rmass;
        vcm.x += velocities[ai].x * rmass;
        vcm.y += velocities[ai].y * rmass;
        vcm.z += velocities[ai].z * rmass;
    }

    vcm.x /= total_mass;
    vcm.y /= total_mass;
    vcm.z /= total_mass;

    for (int ai = 0; ai < ctx.n_atoms; ai++) {
        velocities[ai].x -= vcm.x;
        velocities[ai].y -= vcm.y;
        velocities[ai].z -= vcm.z;
    }
    if (ctx.command_info.requested_gpu) {
        ctx.velocities->upload();
    }
}

void Handler::calc_final_potential(int iteration) {
    ctx.step = iteration;
    reset_energies();
    calc_nonbonded_forces();
    calc_internal_forces(iteration);
    prepare_offdiags(iteration);
}

void Handler::update_energy_totals() {
    if (ctx.command_info.requested_gpu) {
        ctx.energy.download();
        ctx.offdiag_records->download();
    }
    ctx.energy.unpack();
    temperature_->sync_for_output(ctx);  // publish Ukin before combine sums Utot
    ctx.energy.combine(ctx.lambdas->cpu_data_p);
}

void Handler::run() {
    if (!configured_) {
        throw InputError("Handler::run called before configure.");
    }
    if (broken_) {
        throw QgpuError("Previous run failed; call set_state() before running again.");
    }

    try {
        create_outputs();  // moved from initialize(): new files per stage
        init_outputs();

        // 1. temperature calculation
        calc_temperature();

        int num_iterations = ctx.md.steps;
        auto t0 = std::chrono::steady_clock::now();

        for (int i = 0; i < num_iterations; i++) {
            ctx.step = i;
            run_iteration(i);
        }
        auto t1 = std::chrono::steady_clock::now();

        double ms = std::chrono::duration<double, std::milli>(t1 - t0).count();
        printf("MD loop: %.3f ms total, %.4f ms/step (%d steps)\n", ms, ms / num_iterations, num_iterations);

        calc_final_potential(num_iterations);

        finish_outputs();
        shutdown_outputs();
        outputs_.clear();
    } catch (...) {
        broken_ = true;
        outputs_.clear();  // closes files (ofstream destructors)
        throw;
    }
}

void Handler::print_outputs(int iteration) {
    bool needs_energy = false;
    bool needs_restart = false;
    bool needs_trajectory = false;

    for (const auto& output : outputs_) {
        const auto required = output->requirements(ctx, iteration);

        needs_energy |= required.energy;
        needs_restart |= required.restart;
        needs_trajectory |= required.trajectory;
    }

    if (needs_energy) {
        update_energy_totals();
    }

    if (needs_restart) {
        water_boundary_force_->sync_for_output(ctx);
    }

    if (ctx.command_info.requested_gpu) {
        if (needs_trajectory || needs_restart) {
            ctx.coords->download();
        }
        if (needs_restart) {
            ctx.velocities->download();
        }
    }

    for (auto& output : outputs_) {
        output->output(ctx, iteration);
    }
}

void Handler::create_outputs() {
    outputs_.clear();
    outputs_.push_back(std::make_unique<StdOutput>());
    outputs_.push_back(std::make_unique<NativeOutput>(ctx.native_output));
}

void Handler::init_outputs() {
    for (auto& output : outputs_) {
        output->init(ctx);
    }
}

void Handler::finish_outputs() {
    update_energy_totals();
    water_boundary_force_->sync_for_output(ctx);
    if (ctx.command_info.requested_gpu) {
        ctx.coords->download();
        ctx.velocities->download();
    }
    for (auto& output : outputs_) {
        output->finish(ctx);
    }
}

void Handler::shutdown_outputs() {
    for (auto& output : outputs_) {
        output->shutdown();
    }
}

void Handler::reset_energies() {
    ctx.energy.reset();
    ctx.dvelocities->zero();
}

void Handler::calc_internal_forces(int iteration) {
    bonded_force_->calc(ctx);
    restraint_force_->calc(ctx);
    water_boundary_force_->calc(ctx, iteration);
}

void Handler::calc_nonbonded_forces() {
    nonbonded_force_->calc(ctx);
}

void Handler::calc_temperature() {
    temperature_->calc(ctx);
}

void Handler::calc_leapfrog() {
    integrator_->step(ctx);
}

void Handler::prepare_offdiags(int iteration) {
    // offdiags is too easy, so doesn't declare a class to handle it.
    if (ctx.n_offdiags() == 0 || ctx.md.energy <= 0 || iteration <= 0 || iteration % ctx.md.energy != 0) {
        return;
    }

    if (ctx.command_info.requested_gpu) {
        calc_offdiags_cuda(ctx);
    } else {
        calc_offdiags_cpu(ctx);
    }
}
