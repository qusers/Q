#pragma once
#include "bonded_force.h"

class CpuBondedForce : public BondedForce {
   public:
    void calc(Context& ctx) override;

   protected:
    void init_backend(Context& ctx) override {}

   private:
    void calc_bonds(Context& ctx);
    void calc_angles(Context& ctx);
    void calc_torsions(Context& ctx);
    void calc_impropers(Context& ctx);
};