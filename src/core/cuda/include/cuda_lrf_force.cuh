#pragma once
#include "nonbonded_force.h"

struct NonbondedCoordinates {
    const real_t* x;
    const real_t* y;
    const real_t* z;
};

enum GroupPairMode : uint8_t {
    GROUP_PAIR_IGNORE = 0,
    GROUP_PAIR_EXACT = 1,
    GROUP_PAIR_LRF = 2
};

struct ExactEntry {
    int x_start;
    int y_start;
};

// Owns the cutoff/LRF path, including its exact near-field tiles.
// The caller owns the combined atom data and refreshes coordinates before calc().
class CudaLrfForce final {
   public:
    void init(Context& ctx, const NonbondedData& data);
    void calc(Context& ctx, const NonbondedCoordinates& coords);

   private:
    void init_calculation_groups(Context& ctx, const NonbondedCoordinates& coords);
    void init_lrf_coefficients(Context& ctx, const NonbondedCoordinates& coords);
    void calc_exact_tiles(Context& ctx, const NonbondedCoordinates& coords);
    void calc_lrf(Context& ctx, const NonbondedCoordinates& coords);
    void build_exact_atom_tiles(Context& ctx);

    const NonbondedData* data_ = nullptr;

    std::unique_ptr<HostDeviceBuffer<uint8_t>> group_pair_modes_;

    std::unique_ptr<HostDeviceBuffer<ExactEntry>> exact_tiles_;

    std::unique_ptr<HostDeviceBuffer<int>> exact_tile_count_;
    std::unique_ptr<HostDeviceBuffer<int>> list_overflow_;

    std::unique_ptr<HostDeviceBuffer<LrfCoefficients>> lrf_coefficients_;

    size_t exact_tile_capacity_ = 0;

    int n_exact_tiles_ = 0;

    std::unique_ptr<HostDeviceBuffer<int>> lrf_slot_to_group_range_;

    std::unique_ptr<HostDeviceBuffer<uint32_t>> exact_pair_masks_;
};
