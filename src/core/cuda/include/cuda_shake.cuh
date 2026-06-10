#pragma once
#include "host_device_buffer.h"
#include "shake.h"

struct ShakeFastWater {
    int o;
    int h1;
    int h2;
    real_t ra;
    real_t ra_inv;
    real_t rb;
    real_t rc;
    real_t rhh;
    real_t rhh2;
    real_t wo_div_wohh;
    real_t wh_div_wohh;
};

struct ShakeNetwork {
    int center;
    int n_hydrogens;
    int hydrogens[3];  // at most 3
    real_t dist2[3];
    real_t center_winv;
    real_t hydrogen_winv[3];
};

class CudaShake final : public Shake {
   public:
    void apply(Context& ctx) override;
    void initial_shake(Context& ctx) override;
    void cleanup() override;

   protected:
    void init_backend(Context& ctx) override;

   private:
    void apply_to(Context& ctx, coord_t* d_coords, coord_t* d_xcoords);

    bool is_init_backend = false;

    std::unique_ptr<HostDeviceBuffer<ShakeFastWater>> shake_fast_waters;
    std::unique_ptr<HostDeviceBuffer<ShakeNetwork>> shake_networks;

    std::unique_ptr<HostDeviceBuffer<ShakeBond>> fallback_shake_bonds;
    std::vector<int> fallback_color_offsets;
    std::unique_ptr<HostDeviceBuffer<int>> fallback_unconverged;

    void find_shake_fast_water(Context& ctx, std::vector<bool>& optimized);
    void find_shake_network(Context& ctx, std::vector<bool>& optimized);
    void find_fallback_shake_bond(Context& ctx, std::vector<bool>& optimized);
};