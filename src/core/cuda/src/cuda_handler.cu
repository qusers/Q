#include <iostream>

#include "context.h"
#include "cuda/include/cuda_polx_water_force.cuh"
#include "cuda/include/cuda_radix_water_force.cuh"
#include "cuda_bonded_force.cuh"
#include "cuda_handler.cuh"
#include "cuda_integrator.cuh"
#include "cuda_nonbonded_force.cuh"
#include "cuda_restraint_force.cuh"
#include "cuda_shake.cuh"
#include "cuda_shake_v2.cuh"
#include "cuda_temperature.cuh"
#include "init.h"

void CudaHandler::initialize_backend() {
    if (!initialized_) {
        auto parse_result = ctx.init_data_from_files();
        shake_ = create_shake_backend();
        nonbonded_force_ = create_nonbonded_force_backend();
        bonded_force_ = create_bonded_force_backend();
        restraint_force_ = create_restraint_force_backend();
        temperature_ = create_temperature_backend();
        integrator_ = create_integrator_backend();

        ctx.preprocess_data(*shake_);
        nonbonded_force_->init(ctx);
        bonded_force_->init(ctx);
        restraint_force_->init(ctx, parse_result.restrspos, parse_result.restrseqs, parse_result.restrdists, parse_result.restrangs, parse_result.restrwalls);
        temperature_->init(ctx, *shake_);
        integrator_->init(ctx, *shake_, *temperature_);

        init_polx_water_force_kernel_data();
        init_radix_water_force_kernel_data();

        initialized_ = true;
    }
}

void CudaHandler::shutdown() {
    if (initialized_) {
        cleanup_polx_water_force();
        cleanup_radix_water_force();
        initialized_ = false;
    }
}

void CudaHandler::calc_internal_forces(int iteration) {
    auto& host = Context::instance();

    bonded_force_->calc(ctx);
    restraint_force_->calc(ctx);

    if (host.n_waters > 0) {
        calc_radix_water_forces_host(temperature_->data().results->gpu_data_p);
        if (host.md.polarisation) {
            calc_polx_water_forces_host(iteration);
        }
    }
}

void CudaHandler::calc_nonbonded_forces() {
    nonbonded_force_->calc(ctx);
}

void CudaHandler::calc_temperature() {
    temperature_->calc(ctx);
}

void CudaHandler::calc_leapfrog() {
    integrator_->step(ctx);
}

void CudaHandler::reset_energies() {
    Handler::reset_energies();
    Context::instance().cuda_reset_energies();
}

std::unique_ptr<Shake> CudaHandler::create_shake_backend() {
    // return std::make_unique<CudaShakeV2>();
    return std::make_unique<CudaShake>();
}

std::unique_ptr<NonbondedForce> CudaHandler::create_nonbonded_force_backend() {
    return std::make_unique<CudaNonbondedForce>();
}

std::unique_ptr<BondedForce> CudaHandler::create_bonded_force_backend() {
    return std::make_unique<CudaBondedForce>();
}

std::unique_ptr<RestraintForce> CudaHandler::create_restraint_force_backend() {
    return std::make_unique<CudaRestraintForce>();
}

std::unique_ptr<Temperature> CudaHandler::create_temperature_backend() {
    return std::make_unique<CudaTemperature>();
}

std::unique_ptr<Integrator> CudaHandler::create_integrator_backend() {
    return std::make_unique<CudaIntegrator>();
}