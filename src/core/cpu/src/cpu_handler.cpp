#include "cpu_handler.h"

#include "cpu_bonded_force.h"
#include "cpu_integrator.h"
#include "cpu_nonbonded_force.h"
#include "cpu_restraint_force.h"
#include "cpu_shake.h"
#include "cpu_temperature.h"
#include "cpu_water_boundary_force.h"

std::unique_ptr<Shake> CpuHandler::create_shake_backend() {
    return std::make_unique<CpuShake>();
}

std::unique_ptr<NonbondedForce> CpuHandler::create_nonbonded_force_backend() {
    return std::make_unique<CpuNonbondedForce>();
}

std::unique_ptr<BondedForce> CpuHandler::create_bonded_force_backend() {
    return std::make_unique<CpuBondedForce>();
}

std::unique_ptr<RestraintForce> CpuHandler::create_restraint_force_backend() {
    return std::make_unique<CpuRestraintForce>();
}

std::unique_ptr<Temperature> CpuHandler::create_temperature_backend() {
    return std::make_unique<CpuTemperature>();
}

std::unique_ptr<Integrator> CpuHandler::create_integrator_backend() {
    return std::make_unique<CpuIntegrator>();
}
std::unique_ptr<WaterBoundaryForce> CpuHandler::create_water_boundary_force_backend() {
    return std::make_unique<CpuWaterBoundaryForce>();
}