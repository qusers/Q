#pragma once
#include "constraint_force.h"
class CpuShake : public ConstraintForce {
   public:
    void apply(Context& ctx, HostDeviceBuffer<coord_t>& xcoords) override;
    void initial_constraint(Context& ctx) override;

    void init_from_bonds(Context& ctx, const std::vector<ConstraintBond>& bonds) override;

   private:
    void apply_to(Context& ctx, coord_t* coords, coord_t* xcoords);
};