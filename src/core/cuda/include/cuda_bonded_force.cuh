#pragma once
#include "bonded_force.h"

class CudaBondedForce : public BondedForce {
   public:
    void calc(Context& ctx) override;

   protected:
   protected:
    void init_backend(Context& ctx) override; 
};