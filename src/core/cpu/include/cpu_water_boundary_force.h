#pragma once
#include "water_boundary_force.h"

class CpuWaterBoundaryForce : public WaterBoundaryForce {
   public:
    void calc(Context& ctx, int iteration) override;

   protected:
    void init_backend(Context& ctx) override {}  // base init() already allocated the host buffers

   private:
    void calc_radix(Context& ctx);                // was calc_radix_w_forces
    void calc_polx(Context& ctx, int iteration);  // was calc_polx_w_forces
};