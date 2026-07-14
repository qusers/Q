#pragma once
#include "water_boundary_force.h"
class CudaWaterBoundaryForce : public WaterBoundaryForce {
public:
    void calc(Context& ctx, int iteration) override;
    void sync_for_output(Context& ctx) override;

private:
    void calc_radix(Context& ctx);
    void calc_polx(Context& ctx, int iteration);
};