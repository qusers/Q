#pragma once
#include "restraint_force.h"

class CudaRestraintForce : public RestraintForce {
   public:
    void calc(Context& ctx) override;

   protected:
    void init_backend(Context& ctx) override;
};