#pragma once
#include "restraint_force.h"

class CpuRestraintForce : public RestraintForce {
   public:
    void calc(Context& ctx) override;

   protected:
    void init_backend(Context& ctx) override {}

   private:
    void calc_posrestr(Context& ctx);   // pshell + fixed atoms
    void calc_restrpos(Context& ctx);   // atom restraints (FEP)
    void calc_restrseq(Context& ctx);   // sequence restraints
    void calc_restrdis(Context& ctx);   // distance restraints (FEP)
    void calc_restrang(Context& ctx);   // angle restraints (FEP)
    void calc_restrwall(Context& ctx);  // wall restraints
};
