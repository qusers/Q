#include "cpu_nonbonded_qw_force.h"

#include <math.h>

#include "constants.h"
#include "context.h"
#include "vdw_rules.h"
void calc_nonbonded_qw_forces() {
    auto& ctx = Context::instance();
    int i;
    coord_t dO, dH1, dH2;
    double r2O, rH1, rH2, r6O, rO, r2H1, r2H2;
    double dvO, dvH1, dvH2;
    double V_a, V_b, VelO, VelH1, VelH2;
    q_catype_t qi_type;
    double ai_aii, ai_bii;

    if (ctx.A_O == 0) {
        catype_t catype_ow;  // Atom type of first O, H atom

        catype_ow = ctx.catypes[ctx.atypes[ctx.n_atoms_solute].code - 1];

        ctx.A_O = catype_ow.aii_normal;
        ctx.B_O = catype_ow.bii_normal;
    }

    // Loop over O-atoms, q-atoms
    for (int j = ctx.n_atoms_solute; j < ctx.n_atoms; j += 3) {
        for (int qi = 0; qi < ctx.n_qatoms; qi++) {
            i = ctx.q_atoms[qi].a - 1;
            if (ctx.excluded[i] || ctx.excluded[j]) continue;
            dO.x = ctx.coords[j].x - ctx.coords[i].x;
            dO.y = ctx.coords[j].y - ctx.coords[i].y;
            dO.z = ctx.coords[j].z - ctx.coords[i].z;
            dH1.x = ctx.coords[j + 1].x - ctx.coords[i].x;
            dH1.y = ctx.coords[j + 1].y - ctx.coords[i].y;
            dH1.z = ctx.coords[j + 1].z - ctx.coords[i].z;
            dH2.x = ctx.coords[j + 2].x - ctx.coords[i].x;
            dH2.y = ctx.coords[j + 2].y - ctx.coords[i].y;
            dH2.z = ctx.coords[j + 2].z - ctx.coords[i].z;
            r2O = pow(dO.x, 2) + pow(dO.y, 2) + pow(dO.z, 2);
            rH1 = sqrt(1.0 / (pow(dH1.x, 2) + pow(dH1.y, 2) + pow(dH1.z, 2)));
            rH2 = sqrt(1.0 / (pow(dH2.x, 2) + pow(dH2.y, 2) + pow(dH2.z, 2)));
            r6O = r2O * r2O * r2O;
            r2O = 1.0 / r2O;
            rO = sqrt(r2O);
            double r6Oinv = r2O * r2O * r2O;  // 1/r^6 for vdW calculation
            r2H1 = rH1 * rH1;
            r2H2 = rH2 * rH2;

            // Reset potential
            dvO = 0;
            dvH1 = 0;
            dvH2 = 0;

            for (int state = 0; state < ctx.n_lambdas; state++) {
                qi_type = ctx.q_catypes[ctx.q_atypes[qi + ctx.n_qatoms * state].code - 1];

                ai_aii = qi_type.Ai;
                ai_bii = qi_type.Bi;

                if (ctx.topo.vdw_rule == VDW_GEOMETRIC) {
                    calc_vdw_geometric(ai_aii, ctx.A_O, ai_bii, ctx.B_O, r6Oinv, &V_a, &V_b);
                } else {
                    calc_vdw_arithmetic(ai_aii, ctx.A_O, ai_bii, ctx.B_O, r6Oinv, &V_a, &V_b);
                }

                VelO = ctx.topo.coulomb_constant * ctx.crg_ow * ctx.q_charges[qi + ctx.n_qatoms * state].q * rO;
                VelH1 = ctx.topo.coulomb_constant * ctx.crg_hw * ctx.q_charges[qi + ctx.n_qatoms * state].q * rH1;
                VelH2 = ctx.topo.coulomb_constant * ctx.crg_hw * ctx.q_charges[qi + ctx.n_qatoms * state].q * rH2;

                // if (state == 0 && qi == 1) printf("j = %d ai__aii = %f A_O = %f B_O = %f V_a = %f V_b = %f r6O = %f\n", j, ai_aii, A_O, B_O, V_a, V_b, r6O);

                dvO += r2O * (-VelO - (12 * V_a - 6 * V_b)) * ctx.lambdas[state];
                dvH1 -= r2H1 * VelH1 * ctx.lambdas[state];
                dvH2 -= r2H2 * VelH2 * ctx.lambdas[state];

                ctx.EQ_nonbond_qw[state].Ucoul += (VelO + VelH1 + VelH2);
                ctx.EQ_nonbond_qw[state].Uvdw += (V_a - V_b);
            }

            // Note r6O is not the usual 1/rO^6, but rather rO^6. be careful!!!

            // Update forces on Q-atom
            ctx.dvelocities[i].x -= (dvO * dO.x + dvH1 * dH1.x + dvH2 * dH2.x);
            ctx.dvelocities[i].y -= (dvO * dO.y + dvH1 * dH1.y + dvH2 * dH2.y);
            ctx.dvelocities[i].z -= (dvO * dO.z + dvH1 * dH1.z + dvH2 * dH2.z);

            // Update forces on water
            ctx.dvelocities[j].x += dvO * dO.x;
            ctx.dvelocities[j].y += dvO * dO.y;
            ctx.dvelocities[j].z += dvO * dO.z;
            ctx.dvelocities[j + 1].x += dvH1 * dH1.x;
            ctx.dvelocities[j + 1].y += dvH1 * dH1.y;
            ctx.dvelocities[j + 1].z += dvH1 * dH1.z;
            ctx.dvelocities[j + 2].x += dvH2 * dH2.x;
            ctx.dvelocities[j + 2].y += dvH2 * dH2.y;
            ctx.dvelocities[j + 2].z += dvH2 * dH2.z;
        }
    }

#ifdef DEBUG
    printf("q-w: Ecoul = %f Evdw = %f\n", ctx.EQ_nonbond_qw[0].Ucoul, ctx.EQ_nonbond_qw[0].Uvdw);
#endif
}