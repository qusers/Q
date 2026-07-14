#include "handler.h"

#include <chrono>

#include "csv_out.h"
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

    const bool output_step = ctx.md.output > 0 && iteration % ctx.md.output == 0;
    const bool energy_step = ctx.md.energy > 0 && iteration % ctx.md.energy == 0;
    const bool restart_step = iteration % 1000 == 0;
    if (output_step || energy_step) {
        update_energy_totals();  // host energy + temperature for stdout and/or energy file
    }
    if (restart_step) {
        water_boundary_force_->sync_for_output(ctx);
    }
    // 8. print output to files
    print_outputs(iteration);
}

void Handler::initialize() {
    initialize_backend();
    create_outputs();
    init_outputs();
}

void Handler::calc_final_potential(int iteration) {
    reset_energies();
    calc_nonbonded_forces();
    calc_internal_forces(iteration);
    update_energy_totals();
    water_boundary_force_->sync_for_output(ctx);
}

void Handler::run(int num_iterations) {
    // 1. temperature calculation
    calc_temperature();
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
    auto& host = Context::instance();
    if (host.command_info.requested_gpu) {
        host.energy.download();
    }
    host.energy.unpack();
    temperature_->sync_for_output(host);  // publish Ukin before combine sums Utot
    host.energy.combine(host.lambdas->cpu_data_p);
}

void Handler::print_outputs(int iteration) {
    for (auto& output : outputs_) {
        output->output(ctx, iteration);
    }
}

void Handler::create_outputs() {
    outputs_.clear();
    outputs_.push_back(std::make_unique<StdOutput>());

    if (ctx.command_info.input_mode == CommandInputMode::Csv) {
        outputs_.push_back(std::make_unique<CsvOutput>(ctx.command_info.csv_dir));
    } else {
        outputs_.push_back(std::make_unique<NativeOutput>(ctx.native_output));
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

Shake& Handler::shake() {
    return *shake_;
}

NonbondedForce& Handler::nonbonded_force() {
    return *nonbonded_force_;
}

BondedForce& Handler::bonded_force() {
    return *bonded_force_;
}

RestraintForce& Handler::restraint_force() {
    return *restraint_force_;
}

Temperature& Handler::temperature() {
    return *temperature_;
}

Integrator& Handler::integrator() {
    return *integrator_;
}
WaterBoundaryForce& Handler::water_boundary_force() {
    return *water_boundary_force_;
}

void Handler::reset_energies() {
    auto& host = Context::instance();
    host.energy.reset();
    auto& dvelocities = host.dvelocities->cpu_data_p;
    for (int i = 0; i < host.n_atoms; i++) {
        auto& dvel = dvelocities[i];
        dvel.x = 0;
        dvel.y = 0;
        dvel.z = 0;
    }
}
