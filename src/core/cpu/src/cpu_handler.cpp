#include "cpu_handler.h"

#include <time.h>

#include "context.h"
#include "cpu_bonded_force.h"
#include "cpu_integrator.h"
#include "cpu_nonbonded_force.h"
#include "cpu_q_angle_force.h"
#include "cpu_q_bond_force.h"
#include "cpu_q_torsion_force.h"
#include "cpu_restraint_force.h"
#include "cpu_shake.h"
#include "cpu_temperature.h"
#include "cpu_water_boundary_force.h"
#include "init.h"

namespace {

void calc_q_bonded_forces(Context& ctx) {
    for (int state = 0; state < ctx.n_lambdas; state++) {
        calc_qangle_forces(state);
        calc_qbond_forces(state);
        // calc_qtorsion_forces(state);
    }
}

}  // namespace

void CpuHandler::initialize_backend() {
    auto parse_result = ctx.init_data_from_files();
    shake_ = create_shake_backend();
    nonbonded_force_ = create_nonbonded_force_backend();
    bonded_force_ = create_bonded_force_backend();
    restraint_force_ = create_restraint_force_backend();
    temperature_ = create_temperature_backend();
    integrator_ = create_integrator_backend();
    water_boundary_force_ = create_water_boundary_force_backend();
    ctx.preprocess_data(*shake_);
    nonbonded_force_->init(ctx);
    bonded_force_->init(ctx);
    restraint_force_->init(ctx, parse_result.restrspos, parse_result.restrseqs, parse_result.restrdists, parse_result.restrangs, parse_result.restrwalls);
    temperature_->init(ctx, *shake_);
    water_boundary_force_->init(ctx, *temperature_, parse_result.restart_theta_corr);

    integrator_->init(ctx, *shake_, *temperature_);
}

void CpuHandler::shutdown() {
}

void CpuHandler::calc_internal_forces(int iteration) {
    bonded_force_->calc(ctx);
    restraint_force_->calc(ctx);
    water_boundary_force_->calc(ctx, iteration);
}

void CpuHandler::calc_nonbonded_forces() {
    nonbonded_force_->calc(ctx);
}

void CpuHandler::calc_temperature() {
    temperature_->calc(ctx);
}

void CpuHandler::calc_leapfrog() {
    integrator_->step(ctx);
}

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