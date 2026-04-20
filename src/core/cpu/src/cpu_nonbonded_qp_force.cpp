#include "cpu_nonbonded_qp_force.h"

#include <math.h>

#include "constants.h"
#include "context.h"
#include "vdw_rules.h"

void calc_nonbonded_qp_forces() {
    auto& ctx = Context::instance();
    auto &coords = ctx.coords->cpu_data_p;
    auto &dvelocities = ctx.dvelocities->cpu_data_p;
    int i, j;
    coord_t da;
    double r2, r6, r;
    double ai_aii, aj_aii, ai_bii, aj_bii;
    bool bond23, bond14;
    double scaling, Vel, V_a, V_b, dv;

    for (int qi = 0; qi < ctx.n_qatoms; qi++) {
        for (int pj = 0; pj < ctx.n_patoms; pj++) {
            i = ctx.q_atoms[qi];
            j = ctx.p_atoms[pj];

            bond23 = ctx.LJ_matrix[i * ctx.n_atoms_solute + j] == 3;
            bond14 = ctx.LJ_matrix[i * ctx.n_atoms_solute + j] == 1;

            if (bond23) continue;
            if (ctx.excluded[i] || ctx.excluded[j]) continue;

            scaling = bond14 ? ctx.topo.el14_scale : 1;

            da.x = coords[j].x - coords[i].x;
            da.y = coords[j].y - coords[i].y;
            da.z = coords[j].z - coords[i].z;

            r2 = pow(da.x, 2) + pow(da.y, 2) + pow(da.z, 2);

            r6 = r2 * r2 * r2;
            r2 = 1 / r2;
            r = sqrt(r2);
            double r6inv = r2 * r2 * r2;  // 1/r^6 for vdW calculation

            for (int state = 0; state < ctx.n_lambdas; state++) {
                const catype_t& qi_type = ctx.unified_catype(i, state);
                const catype_t& aj_type = ctx.unified_catype(j, state);

                ai_aii = bond14 ? qi_type.aii_1_4 : qi_type.aii_normal;
                aj_aii = bond14 ? aj_type.aii_1_4 : aj_type.aii_normal;
                ai_bii = bond14 ? qi_type.bii_1_4 : qi_type.bii_normal;
                aj_bii = bond14 ? aj_type.bii_1_4 : aj_type.bii_normal;

                Vel = ctx.topo.coulomb_constant * scaling * ctx.unified_ccharge(i, state).charge * ctx.unified_ccharge(j, state).charge * r;
                if (ctx.topo.vdw_rule == VDW_GEOMETRIC) {
                    calc_vdw_geometric(ai_aii, aj_aii, ai_bii, aj_bii, r6inv, &V_a, &V_b);
                } else {
                    calc_vdw_arithmetic(ai_aii, aj_aii, ai_bii, aj_bii, r6inv, &V_a, &V_b);
                }
                dv = r2 * (-Vel - (12 * V_a - 6 * V_b)) * ctx.lambdas[state];

                // Update forces
                dvelocities[i].x -= dv * da.x;
                dvelocities[i].y -= dv * da.y;
                dvelocities[i].z -= dv * da.z;
                dvelocities[j].x += dv * da.x;
                dvelocities[j].y += dv * da.y;
                dvelocities[j].z += dv * da.z;

                // Update Q totals
                ctx.EQ_nonbond_qp[state].Ucoul += Vel;
                ctx.EQ_nonbond_qp[state].Uvdw += (V_a - V_b);
            }
        }
    }
}
