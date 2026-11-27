#pragma once
#include "integrator.h"
class CpuIntegrator : public Integrator {
   public:
    void step(Context& ctx) override;  // body = today's cpu calc_leapfrog
};