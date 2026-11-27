#pragma once

#include "temperature.h"

class CpuTemperature : public Temperature {
   public:
    void calc(Context& ctx) override;
};