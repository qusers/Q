#pragma once

#include "handler.h"

class CpuHandler : public Handler {
   protected:
    std::unique_ptr<Shake> create_shake_backend() override;
    std::unique_ptr<NonbondedForce> create_nonbonded_force_backend() override;
    std::unique_ptr<BondedForce> create_bonded_force_backend() override;
    std::unique_ptr<RestraintForce> create_restraint_force_backend() override;
    std::unique_ptr<Temperature> create_temperature_backend() override;
    std::unique_ptr<Integrator> create_integrator_backend() override;
    std::unique_ptr<WaterBoundaryForce> create_water_boundary_force_backend() override;
};
