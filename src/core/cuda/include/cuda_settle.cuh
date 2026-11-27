#pragma once
#include <cstdint>
#include <memory>

#include "constraint_force.h"
#include "host_device_buffer.h"
#include "md_types.h"
struct SettleFastWater {
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

SettleFastWater make_fast_water(int o, int h1, int h2, double roh2, double rhh2, double wo, double wh);

class CudaSettleSolver {
   public:
    void init(const std::vector<SettleFastWater>& settle_fast_waters);
    void apply(coord_t* d_coords, const coord_t* d_xcoords, const double* d_winv);

   private:
    std::unique_ptr<HostDeviceBuffer<SettleFastWater>> settle_fast_waters_;
};

class CudaSettle final : public ConstraintForce {
   public:
    void apply(Context& ctx, HostDeviceBuffer<coord_t>& xcoords) override;
    void initial_constraint(Context& ctx) override;
    void init_from_bonds(Context& ctx, const std::vector<ConstraintBond>& bonds) override;

   protected:
    void init_backend(Context& ctx) override;

   private:
    void apply_to(Context& ctx, coord_t* d_coords, coord_t* d_xcoords);

    CudaSettleSolver solver_;
};
