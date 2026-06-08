#pragma once
#include "shake.h"

struct ShakeFastWater {
    int o;
    int h1;
    int h2;
    real_t ra;
    real_t ra_inv;
    real_t rb;
    real_t rc;
    real_t rc2;
    real_t hhhh;
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
};