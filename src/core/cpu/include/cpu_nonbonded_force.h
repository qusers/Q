#pragma once
#include "nonbonded_force.h"

class CpuNonbondedForce final : public NonbondedForce {
   public:
    void calc(Context& ctx) override;

};