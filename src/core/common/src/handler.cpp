#include "handler.h"

#include <chrono>
#include <unordered_set>

#include "native_out.h"
#include "std_output.h"

void Handler::run_iteration(int iteration) {
    reset_energies();

    // 3. nonbonded forces
    calc_nonbonded_forces();

    // 2. bonded forces and some constraints calculations
    calc_internal_forces(iteration);

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
    std::vector<ParseResult> replica_results = ctx.get_parse_results();
    ctx.init(replica_results);
    const ParseResult& parsed = replica_results.front();
    shake_ = create_shake_backend();
    nonbonded_force_ = create_nonbonded_force_backend();
    bonded_force_ = create_bonded_force_backend();
    restraint_force_ = create_restraint_force_backend();
    temperature_ = create_temperature_backend();
    integrator_ = create_integrator_backend();
    water_boundary_force_ = create_water_boundary_force_backend();
    shake_->init(ctx, parsed);
    bonded_force_->init(ctx, parsed, shake_->data());
    nonbonded_force_->init(ctx);
    restraint_force_->init(ctx, parsed.restrspos, parsed.restrseqs, parsed.restrdists, parsed.restrangs, parsed.restrwalls);
    temperature_->init(ctx, *shake_);
    water_boundary_force_->init(ctx, *temperature_, replica_results);
    integrator_->init(ctx, *shake_, *temperature_);

    if (ctx.fresh_start) {
        shake_->initial_shake(ctx);
        stop_cm_translation();
    }

    create_outputs();
    init_outputs();
}

void Handler::stop_cm_translation() {
    auto& atypes = ctx.atypes->cpu_data_p;
    auto& catypes = ctx.catypes->cpu_data_p;
    auto& velocities = ctx.velocities->cpu_data_p;
    double total_mass = 0;

    for (int ai = 0; ai < ctx.n_atoms; ai++) {
        const double rmass = catypes[atypes[ai].code - 1].m;
        total_mass += rmass;
    }

    for (int replica = 0; replica < ctx.n_replicates(); replica++) {
        const int atom_offset = replica * ctx.n_atoms;
        coord_t vcm = {};
        for (int ai = 0; ai < ctx.n_atoms; ai++) {
            const int atom = atom_offset + ai;
            const double mass = catypes[atypes[ai].code - 1].m;
            vcm.x += velocities[atom].x * mass;
            vcm.y += velocities[atom].y * mass;
            vcm.z += velocities[atom].z * mass;
        }
        vcm.x /= total_mass;
        vcm.y /= total_mass;
        vcm.z /= total_mass;

        for (int ai = 0; ai < ctx.n_atoms; ai++) {
            const int atom = atom_offset + ai;

            velocities[atom].x -= vcm.x;
            velocities[atom].y -= vcm.y;
            velocities[atom].z -= vcm.z;
        }
    }

    if (ctx.command_info.requested_gpu) {
        ctx.velocities->upload();
    }
}

void Handler::calc_final_potential(int iteration) {
    reset_energies();
    calc_nonbonded_forces();
    calc_internal_forces(iteration);
    update_energy_totals();
    water_boundary_force_->sync_for_output(ctx);
}

void Handler::run() {
    // 1. temperature calculation
    calc_temperature();

    int num_iterations = ctx.md.steps;
    auto t0 = std::chrono::steady_clock::now();

    for (int i = 0; i < num_iterations; i++) {
        run_iteration(i);
    }
    auto t1 = std::chrono::steady_clock::now();

    double ms = std::chrono::duration<double, std::milli>(t1 - t0).count();
    printf("MD loop: %.3f ms total, %.4f ms/step (%d steps)\n", ms, ms / num_iterations, num_iterations);

    calc_final_potential(num_iterations);

    finish_outputs();
    shutdown_outputs();
    outputs_.clear();
}

void Handler::update_energy_totals() {
    if (ctx.command_info.requested_gpu) {
        ctx.energy.download();
    }
    for (int replica = 0; replica < ctx.n_replicates(); replica++) {
        ctx.energy.unpack(replica);
        temperature_->sync_for_output(ctx, replica);
        ctx.energy.combine(ctx.lambdas->cpu_data_p, replica);
    }
}

void Handler::print_outputs(int iteration) {
    bool needs_energy = false;
    bool needs_restart = false;

    for (const auto& output : outputs_) {
        const auto required = output->requirements(ctx, iteration);

        needs_energy |= required.energy;
        needs_restart |= required.restart;
    }

    if (needs_energy) {
        update_energy_totals();
    }

    if (needs_restart) {
        water_boundary_force_->sync_for_output(ctx);
    }

    for (auto& output : outputs_) {
        output->output(ctx, iteration);
    }
}

void Handler::create_outputs() {
    outputs_.clear();
    outputs_.push_back(std::make_unique<StdOutput>());

    std::unordered_set<std::string> output_files;
    for (int replica = 0; replica < ctx.n_replicates(); replica++) {
        NativeOutputConfig config = ctx.native_outputs.at(replica);

        auto already_used = [&](const std::string& path) {
            return !path.empty() && output_files.count(path) > 0;
        };
        const bool collision = already_used(config.final_file) || already_used(config.trajectory_file) || already_used(config.energy_file);

        if (collision) {
            std::fprintf(stderr, ">>> WARNING: native output for replica %d is skipped because an output filename is already used by an earlier replica.\n", replica);
            continue;
        }

        auto remember = [&](const std::string& path) {
            if (!path.empty()) {
                output_files.insert(path);
            }
        };

        remember(config.final_file);
        remember(config.trajectory_file);
        remember(config.energy_file);

        outputs_.push_back(std::make_unique<NativeOutput>(std::move(config), replica));
    }
}

void Handler::init_outputs() {
    for (auto& output : outputs_) {
        output->init(ctx);
    }
}

void Handler::finish_outputs() {
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
