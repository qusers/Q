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
    int x_len;  // <= 32
    int y_start;
    int y_len;  // <= 32
    uint8_t diagonal;
};

struct LrfPairEntry {
    int range1;
    int range2;

    int group1;
    int group2;
};

class CudaNonbondedForce final : public NonbondedForce {
   public:
    void calc(Context& ctx) override;

   protected:
    void init_backend(Context& ctx) override;

   private:
    void calc_all_direct_pairs(Context& ctx);
    void init_calculation_groups(Context& ctx);
    void init_calculation_groups_by_switch(Context& ctx);
    void init_calculation_groups_by_all_atoms(Context& ctx);
    void init_lrf_coefficients(Context &ctx);
    void calc_exact_tiles(Context& ctx);
    void calc_lrf(Context &ctx);

    std::unique_ptr<HostDeviceBuffer<real_t>> coord_x_, coord_y_, coord_z_;

    std::unique_ptr<HostDeviceBuffer<uint8_t>> group_pair_modes_;

    std::unique_ptr<HostDeviceBuffer<ExactEntry>> exact_tiles_;
    std::unique_ptr<HostDeviceBuffer<LrfPairEntry>> lrf_group_pairs_;

    std::unique_ptr<HostDeviceBuffer<int>> exact_tile_count_;
    std::unique_ptr<HostDeviceBuffer<int>> lrf_pair_count_;
    std::unique_ptr<HostDeviceBuffer<int>> list_overflow_;

    std::unique_ptr<HostDeviceBuffer<LrfCoefficients>> lrf_coefficients_;

    size_t exact_tile_capacity_ = 0;
    size_t lrf_pair_capacity_ = 0;

    int n_exact_tiles_ = 0;
    int n_lrf_pairs_ = 0;
};
