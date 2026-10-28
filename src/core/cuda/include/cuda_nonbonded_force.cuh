#pragma once
#include "cuda_lrf_force.cuh"
#include "nonbonded_force.h"

class CudaNonbondedForce final : public NonbondedForce {
   public:
    void calc(Context& ctx) override;

   protected:
    void init_backend(Context& ctx) override;

   private:
    void calc_all_direct_pairs(Context& ctx);

    std::unique_ptr<HostDeviceBuffer<real_t>> coord_x_, coord_y_, coord_z_;
    std::unique_ptr<CudaLrfForce> lrf_;
};
