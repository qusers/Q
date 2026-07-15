#pragma once
#include "shake.h"
class CpuShake : public Shake {
   public:
    void apply(Context& ctx, HostDeviceBuffer<coord_t>& xcoords) override;
    void initial_shake(Context& ctx, HostDeviceBuffer<coord_t>& xcoords_buffer) override;

   private:
    void apply_to(Context& ctx, coord_t* coords, coord_t* xcoords);
};