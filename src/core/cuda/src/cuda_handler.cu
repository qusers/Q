#include <iostream>

#include "cuda_bonded_force.cuh"
#include "cuda_handler.cuh"
#include "cuda_integrator.cuh"
#include "cuda_nonbonded_force.cuh"
#include "cuda_restraint_force.cuh"
#include "cuda_shake.cuh"
#include "cuda_shake_v2.cuh"
#include "cuda_temperature.cuh"
#include "cuda_water_boundary_force.cuh"

std::unique_ptr<Shake> CudaHandler::create_shake_backend() {
    // return std::make_unique<CudaShakeV2>();
    return std::make_unique<CudaShake>(true);
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
std::unique_ptr<WaterBoundaryForce> CudaHandler::create_water_boundary_force_backend() {
    return std::make_unique<CudaWaterBoundaryForce>();
}
