#pragma once
#include "nonbonded_force.h"
#include "precision.h"

class CudaNonbondedForce final : public NonbondedForce {
   public:
    void calc(Context& ctx) override;

   protected:
    void init_backend(Context& ctx) override;

   private:
    std::unique_ptr<HostDeviceBuffer<energy_accum_t>> e_coul_;
    std::unique_ptr<HostDeviceBuffer<energy_accum_t>> e_vdw_;
    std::unique_ptr<HostDeviceBuffer<real_t3>> coords;
};
