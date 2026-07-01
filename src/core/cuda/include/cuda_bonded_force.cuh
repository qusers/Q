#pragma once
#include "bonded_force.h"

class CudaBondedForce : public BondedForce {
   public:
    void calc(Context& ctx) override;

   protected:
   protected:
    void init_backend(Context& ctx) override; 

   private:
    // one accumulator buffer per kind, length bonded_region_slots() = 2 (P, W)
    std::unique_ptr<HostDeviceBuffer<energy_accum_t>> e_bond_, e_angle_, e_torsion_, e_improper_;
};