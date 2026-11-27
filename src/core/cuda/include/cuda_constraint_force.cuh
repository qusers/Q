#pragma once
#include "constraint_force.h"

class CudaConstraintForce final : public ConstraintForce {
    void apply(Context& ctx, HostDeviceBuffer<coord_t>& xcoords) override;
    void initial_constraint(Context& ctx) override;
    void cleanup() override;
    void init_from_bonds(Context& ctx, const std::vector<ConstraintBond>& bonds) override;

   protected:
    void init_backend(Context& ctx) override;

    std::unique_ptr<ConstraintForce> solute_constraint_force_;
    std::unique_ptr<ConstraintForce> solvent_constraint_force_;
    std::unique_ptr<ConstraintForce> common_constraint_force_;
};
