#include "cpu_nonbonded_qw_force.h"
#include "cpu_force_accumulation.h"

#include <math.h>
#include <vector>

#include "constants.h"
#include "context.h"
#include "vdw_rules.h"
void calc_nonbonded_qw_forces() {
    auto& ctx = Context::instance();
    auto *lambdas = ctx.lambdas->cpu_data_p;
    auto &coords = ctx.coords->cpu_data_p;
    auto &dvelocities = ctx.dvelocities->cpu_data_p;
    auto *excluded = ctx.excluded->cpu_data_p;
    int i;
    coord_t dO, dH1, dH2;
    double r2O, rH1, rH2, rO, r2H1, r2H2;
    double dvO, dvH1, dvH2;
    double V_a, V_b, VelO, VelH1, VelH2;
    double ai_aii, ai_bii;
    std::vector<energy_accum_t> ucoul(ctx.n_lambdas, 0);
    std::vector<energy_accum_t> uvdw(ctx.n_lambdas, 0);

    // Loop over O-atoms, q-atoms
    for (int j = ctx.n_atoms_solute; j < ctx.n_atoms; j += 3) {
        const catype_t& ow_type = ctx.unified_catype(j, 0);
        double ow_charge = ctx.unified_ccharge(j, 0).charge;
        double hw1_charge = ctx.unified_ccharge(j + 1, 0).charge;
        double hw2_charge = ctx.unified_ccharge(j + 2, 0).charge;
        for (int qi = 0; qi < ctx.n_qatoms; qi++) {
            i = ctx.q_atoms[qi];
            if (excluded[i] || excluded[j]) continue;
            dO.x = coords[j].x - coords[i].x;
            dO.y = coords[j].y - coords[i].y;
            dO.z = coords[j].z - coords[i].z;
            dH1.x = coords[j + 1].x - coords[i].x;
            dH1.y = coords[j + 1].y - coords[i].y;
            dH1.z = coords[j + 1].z - coords[i].z;
            dH2.x = coords[j + 2].x - coords[i].x;
            dH2.y = coords[j + 2].y - coords[i].y;
            dH2.z = coords[j + 2].z - coords[i].z;
            r2O = dO.x * dO.x + dO.y * dO.y + dO.z * dO.z;
            rH1 = std::sqrt(1.0 / (dH1.x * dH1.x + dH1.y * dH1.y + dH1.z * dH1.z));
            rH2 = std::sqrt(1.0 / (dH2.x * dH2.x + dH2.y * dH2.y + dH2.z * dH2.z));
            r2O = 1.0 / r2O;
            rO = std::sqrt(r2O);
            const double r6Oinv = r2O * r2O * r2O;  // 1/r^6 for vdW calculation
            r2H1 = rH1 * rH1;
            r2H2 = rH2 * rH2;

            // Reset potential
            dvO = 0;
            dvH1 = 0;
            dvH2 = 0;

            for (int state = 0; state < ctx.n_lambdas; state++) {
                const catype_t& qi_type = ctx.unified_catype(i, state);

                ai_aii = qi_type.aii_normal;
                ai_bii = qi_type.bii_normal;

                if (ctx.topo.vdw_rule == VDW_GEOMETRIC) {
                    calc_vdw_geometric(ai_aii, ow_type.aii_normal, ai_bii, ow_type.bii_normal, r6Oinv, &V_a, &V_b);
                } else {
                    calc_vdw_arithmetic(ai_aii, ow_type.aii_normal, ai_bii, ow_type.bii_normal, r6Oinv, &V_a, &V_b);
                }

                const double q_charge = ctx.unified_ccharge(i, state).charge;
                VelO =  ow_charge * q_charge * rO * ctx.topo.coulomb_constant;
                VelH1 = hw1_charge * q_charge * rH1 * ctx.topo.coulomb_constant;
                VelH2 = hw2_charge * q_charge * rH2 * ctx.topo.coulomb_constant;


                const double lambda = lambdas[state];
                dvO += r2O * (-VelO - (12.0 * V_a - 6.0 * V_b)) * lambda;
                dvH1 -= r2H1 * VelH1 * lambda;
                dvH2 -= r2H2 * VelH2 * lambda;

                add_energy(ucoul[state], VelO + VelH1 + VelH2);
                add_energy(uvdw[state], V_a - V_b);
            }

            // Note r6O is not the usual 1/rO^6, but rather rO^6. be careful!!!

            // Update forces on Q-atom
            add_force(dvelocities[i].x, -(dvO * dO.x + dvH1 * dH1.x + dvH2 * dH2.x));
            add_force(dvelocities[i].y, -(dvO * dO.y + dvH1 * dH1.y + dvH2 * dH2.y));
            add_force(dvelocities[i].z, -(dvO * dO.z + dvH1 * dH1.z + dvH2 * dH2.z));

            // Update forces on water
            add_force(dvelocities[j].x, dvO * dO.x);
            add_force(dvelocities[j].y, dvO * dO.y);
            add_force(dvelocities[j].z, dvO * dO.z);


            add_force(dvelocities[j + 1].x, dvH1 * dH1.x);
            add_force(dvelocities[j + 1].y, dvH1 * dH1.y);
            add_force(dvelocities[j + 1].z, dvH1 * dH1.z);

            add_force(dvelocities[j + 2].x, dvH2 * dH2.x);
            add_force(dvelocities[j + 2].y, dvH2 * dH2.y);
            add_force(dvelocities[j + 2].z, dvH2 * dH2.z);
        }
    }

    for (int state = 0; state < ctx.n_lambdas; state++) {
        ctx.EQ_nonbond_qw[state].Ucoul += energy_from_accum(ucoul[state]);
        ctx.EQ_nonbond_qw[state].Uvdw += energy_from_accum(uvdw[state]);
    }
}
