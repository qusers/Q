#pragma once
#include "constraint_force.h"
#include "lincs.h"

struct LincsData {
    std::vector<double> bond_inv_sqrt_inv_mass_sum;
    std::vector<std::vector<int>> bond_graph;
    std::vector<std::vector<double>> normalized_c;
};

class CpuLincs : public ConstraintForce {
   public:
    explicit CpuLincs(LincsSettings lincs_settings = {}) : lincs_settings_(lincs_settings) {
    }

    void apply(Context& ctx, HostDeviceBuffer<coord_t>& xcoords) override;
    void initial_constraint(Context& ctx) override;
    void init_from_bonds(Context& ctx, const std::vector<ConstraintBond>& bonds) override;

   protected:
    void init_backend(Context& ctx) override;

   private:
    void apply_to(Context& ctx, coord_t* coords, const coord_t* xcoords);
    LincsData lincs_data_;
    LincsSettings lincs_settings_;

    int get_sign(const int a, const ConstraintBond& bond);

    void solve(Context& ctx, std::vector<double> errors, coord_t* coords, const coord_t* xcoords);
    bool check_ready(const coord_t* coords);
};