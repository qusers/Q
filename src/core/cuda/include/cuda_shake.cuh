#pragma once
#include "host_device_buffer.h"
#include "shake.h"

struct ShakeFastWater {
    int o;
    int h1;
    int h2;
    double ra;
    double ra_inv;
    double rb;
    double rc;
    double rhh;
    double rhh2;
    double wo_div_wohh;
    double wh_div_wohh;
};

struct ShakeNetwork {
    int center;
    int n_hydrogens;
    int hydrogens[3];  // at most 3
    double dist2[3];
    double center_winv;
    double hydrogen_winv[3];
};

class CudaShake final : public Shake {
   public:
    void apply(Context& ctx, HostDeviceBuffer<coord_t>& xcoords) override;
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
    std::unique_ptr<HostDeviceBuffer<int>> fallback_color_offsets;
    std::unique_ptr<HostDeviceBuffer<int>> fallback_unconverged;
    int fallback_n_colors = 0;
    int fallback_coop_blocks = 0;

    void find_shake_fast_water(Context& ctx, std::vector<bool>& optimized);
    void find_shake_network(Context& ctx, std::vector<bool>& optimized);
    void find_fallback_shake_bond(Context& ctx, std::vector<bool>& optimized);
};