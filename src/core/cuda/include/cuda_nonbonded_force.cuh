#pragma once
#include "nonbonded_force.h"
#include "precision.h"

enum GroupPairMode : uint8_t {
    GROUP_PAIR_IGNORE = 0,
    GROUP_PAIR_EXACT = 1,
    GROUP_PAIR_LRF = 2
};

struct ExactEntry {
    int x_start;
    int x_len;  // <= 128
    int y_start;
    int y_len;  // <= 32
};

class CudaNonbondedForce final : public NonbondedForce {
   public:
    void calc(Context& ctx) override;

   protected:
    void init_backend(Context& ctx) override;

   private:
    std::unique_ptr<HostDeviceBuffer<real_t>> coord_x, coord_y, coord_z;
    void calc_all_direct_pairs(Context& ctx);
    void init_calculation_groups(Context& ctx);
};
