#include <iostream>

#include "common/include/context.h"
#include "cuda/include/cuda_handler.cuh"
#include "cuda/include/cuda_leapfrog.cuh"
#include "cuda/include/cuda_polx_water_force.cuh"
#include "cuda/include/cuda_pshell_force.cuh"
#include "cuda/include/cuda_radix_water_force.cuh"
#include "cuda/include/cuda_temperature.cuh"
#include "cuda_bonded_force.cuh"
#include "cuda_nonbonded_force.cuh"
#include "cuda_restraint_force.cuh"
#include "cuda_shake.cuh"
#include "cuda_shake_v2.cuh"
#include "init.h"

void CudaHandler::initialize_backend() {
    if (!initialized_) {
        shake_ = create_shake_backend();
        nonbonded_force_ = create_nonbonded_force_backend();
        bonded_force_ = create_bonded_force_backend();
        restraint_force_ = create_restraint_force_backend();

        ctx.preprocess_data(*shake_);
        nonbonded_force_->init(ctx);
        bonded_force_->init(ctx);
        restraint_force_->init(ctx);

        init_leapfrog_kernel_data();
        init_polx_water_force_kernel_data();
        init_pshell_force_kernel_data();
        init_radix_water_force_kernel_data();
        init_temperature_kernel_data();

        initialized_ = true;
    }
}

void CudaHandler::shutdown() {
    if (initialized_) {
        cleanup_leapfrog();
        cleanup_polx_water_force();
        cleanup_pshell_force();
        cleanup_radix_water_force();
        cleanup_temperature();
        initialized_ = false;
    }
}

void CudaHandler::calc_internal_forces(int iteration) {
    auto& host = Context::instance();

    bonded_force_->calc(ctx);
    restraint_force_->calc(ctx);

    if (host.n_waters > 0) {
        calc_radix_water_forces_host();
        if (host.md.polarisation) {
            calc_polx_water_forces_host(iteration);
        }
    }
}

void CudaHandler::calc_nonbonded_forces() {
    nonbonded_force_->calc(ctx);
}

void CudaHandler::calc_temperature() {
    ::calc_temperature_host();
}

void CudaHandler::calc_leapfrog() {
    calc_leapfrog_host(*shake_);
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