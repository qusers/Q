#pragma once
#include "temperature.h"

class CudaTemperature : public Temperature {
   public:
    void calc(Context& ctx) override;
    void sync_for_output(Context& ctx) override;

   protected:
    void init_backend(Context& ctx) override;

   private:
    std::unique_ptr<HostDeviceBuffer<energy_accum_t>> accum_;  // 6 reduction slots
};