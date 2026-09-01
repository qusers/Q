#pragma once
#include <vector>

#include "nonbonded_force.h"

class CpuNonbondedForce final : public NonbondedForce {
   public:
    void calc(Context& ctx) override;

   protected:
    void init_backend(Context& ctx) override;

   private:
    void calc_all_direct_pairs(Context& ctx);
    void init_calculation_groups(Context& ctx);
    void init_exact_atom_pairs(Context& ctx);
    void calc_direct_pair(Context& ctx, int slot1, int slot2);
    void calc_exact_pairs(Context& ctx);
    void init_lrf_coefficients(Context& ctx);
    void calc_lrf(Context& ctx);

    std::vector<std::pair<int, int>> exact_calculation_groups_, lrf_calculation_groups_;
    std::vector<std::pair<int, int>> exact_atom_pairs_;
    std::vector<LrfCoefficients> lrf_coefficients_;
    std::vector<int> non_q_slot_by_atom_;
};