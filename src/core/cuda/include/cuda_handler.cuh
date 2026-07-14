#pragma once

#include "handler.h"

class CudaHandler : public Handler {
   public:
    // Release device resources.
    void shutdown() override;

   protected:
    bool initialized_ = false;

    void initialize_backend() override;
    void calc_internal_forces(int iteration) override;
    void calc_nonbonded_forces() override;
    void calc_temperature() override;
    void calc_leapfrog() override;

    void reset_energies() override;
    std::unique_ptr<Shake> create_shake_backend() override;
    std::unique_ptr<NonbondedForce> create_nonbonded_force_backend() override;
    std::unique_ptr<BondedForce> create_bonded_force_backend() override;
    std::unique_ptr<RestraintForce> create_restraint_force_backend() override;
    std::unique_ptr<Temperature> create_temperature_backend() override;
    std::unique_ptr<Integrator> create_integrator_backend() override;
    std::unique_ptr<WaterBoundaryForce> create_water_boundary_force_backend() override;
};
