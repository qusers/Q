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
    int y_start;
};

class CudaNonbondedForce final : public NonbondedForce {
   public:
    void calc(Context& ctx) override;

   protected:
    void init_backend(Context& ctx) override;

   private:
    void calc_all_direct_pairs(Context& ctx);
    void init_calculation_groups(Context& ctx);
    void init_lrf_coefficients(Context& ctx);
    void calc_exact_tiles(Context& ctx);
    void calc_lrf(Context& ctx);
    void build_exact_atom_tiles(Context& ctx);

    std::unique_ptr<HostDeviceBuffer<real_t>> coord_x_, coord_y_, coord_z_;

    std::unique_ptr<HostDeviceBuffer<uint8_t>> group_pair_modes_;

    std::unique_ptr<HostDeviceBuffer<ExactEntry>> exact_tiles_;

    std::unique_ptr<HostDeviceBuffer<int>> exact_tile_count_;
    std::unique_ptr<HostDeviceBuffer<int>> list_overflow_;

    std::unique_ptr<HostDeviceBuffer<LrfCoefficients>> lrf_coefficients_;

    size_t exact_tile_capacity_ = 0;

    int n_exact_tiles_ = 0;

    std::unique_ptr<HostDeviceBuffer<int>> lrf_slot_to_group_range_;

    std::unique_ptr<HostDeviceBuffer<uint32_t>> exact_pair_masks_;
    std::unique_ptr<HostDeviceBuffer<uint32_t>> exact_pair_14_masks_;
};
