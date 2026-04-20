#include "cpu_nonbonded_qq_force.h"

#include <math.h>

#include "constants.h"
#include "context.h"
#include "vdw_rules.h"

void calc_nonbonded_qq_forces() {
    auto& ctx = Context::instance();
    auto &coords = ctx.coords->cpu_data_p;
    auto &dvelocities = ctx.dvelocities->cpu_data_p;
    int ai, aj;
    double crg_i, crg_j;
    double elscale, scaling;
    bool bond23, bond14;
    coord_t da;
    double r2a, ra, r6a;
    double Vela, V_a, V_b;
    double dva;
    double ai_aii, aj_aii, ai_bii, aj_bii;

    for (int state = 0; state < ctx.n_lambdas; state++) {
        for (int qi = 0; qi < ctx.n_qatoms; qi++) {
            for (int qj = qi + 1; qj < ctx.n_qatoms; qj++) {
                ai = ctx.q_atoms[qi];
                aj = ctx.q_atoms[qj];

                crg_i = ctx.unified_ccharge(ai, state).charge;
                crg_j = ctx.unified_ccharge(aj, state).charge;

                bond23 = ctx.LJ_matrix[ai * ctx.n_atoms_solute + aj] == 3;
                bond14 = ctx.LJ_matrix[ai * ctx.n_atoms_solute + aj] == 1;

                if (bond23) continue;
                if (ctx.excluded[ai] || ctx.excluded[aj]) continue;

                scaling = bond14 ? ctx.topo.el14_scale : 1;

                elscale = 1;
                for (int k = 0; k < ctx.n_qelscales; k++) {
                    if (ctx.q_elscales[k + ctx.n_qelscales * state].qi == qi + 1 && ctx.q_elscales[k + ctx.n_qelscales * state].qj == qj + 1) {
                        elscale = ctx.q_elscales[k + ctx.n_qelscales * state].mu;
                    }
                }

                const catype_t& qi_type = ctx.unified_catype(ai, state);
                const catype_t& qj_type = ctx.unified_catype(aj, state);

                da.x = coords[aj].x - coords[ai].x;
                da.y = coords[aj].y - coords[ai].y;
                da.z = coords[aj].z - coords[ai].z;
                r2a = 1 / (pow(da.x, 2) + pow(da.y, 2) + pow(da.z, 2));
                ra = sqrt(r2a);
                r6a = r2a * r2a * r2a;

                Vela = scaling * ctx.topo.coulomb_constant * crg_i * crg_j * ra * elscale;

                ai_aii = bond14 ? qi_type.aii_1_4 : qi_type.aii_normal;
                aj_aii = bond14 ? qj_type.aii_1_4 : qj_type.aii_normal;
                ai_bii = bond14 ? qi_type.bii_1_4 : qi_type.bii_normal;
                aj_bii = bond14 ? qj_type.bii_1_4 : qj_type.bii_normal;

                if (ctx.topo.vdw_rule == VDW_GEOMETRIC) {
                    calc_vdw_geometric(ai_aii, aj_aii, ai_bii, aj_bii, r6a, &V_a, &V_b);
                } else {
                    calc_vdw_arithmetic(ai_aii, aj_aii, ai_bii, aj_bii, r6a, &V_a, &V_b);
                }
                dva = r2a * (-Vela - 12 * V_a + 6 * V_b) * ctx.lambdas[state];

                dvelocities[ai].x -= dva * da.x;
                dvelocities[ai].y -= dva * da.y;
                dvelocities[ai].z -= dva * da.z;

                dvelocities[aj].x += dva * da.x;
                dvelocities[aj].y += dva * da.y;
                dvelocities[aj].z += dva * da.z;

                ctx.EQ_nonbond_qq[state].Ucoul += Vela;
                ctx.EQ_nonbond_qq[state].Uvdw += (V_a - V_b);
            }
        }
    }

#ifdef DEBUG
    printf("q-q: Ecoul = %f Evdw = %f\n", ctx.EQ_nonbond_qq[0].Ucoul, ctx.EQ_nonbond_qq[0].Uvdw);
#endif
}
