#pragma once
#include "shake.h"
class CpuShake final : public Shake {
   public:
    void apply(Context& ctx) override;
    void initial_shake(Context& ctx) override;

   private:
    void apply_to(Context& ctx, coord_t* coords, coord_t* xcoords);
};