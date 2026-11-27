#pragma once
#include "integrator.h"
class CudaIntegrator : public Integrator {
   public:
    void step(Context& ctx) override;

};