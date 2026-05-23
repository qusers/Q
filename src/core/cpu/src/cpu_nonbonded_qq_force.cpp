#include "cpu_nonbonded_qq_force.h"

#include <math.h>

#include "constants.h"
#include "context.h"
#include "coulomb.h"
#include "softcore.h"
#include "vdw_rules.h"

void calc_nonbonded_qq_forces() {
    auto& ctx = Context::instance();
    auto *lambdas = ctx.lambdas->cpu_data_p;
    auto &coords = ctx.coords->cpu_data_p;
    auto &dvelocities = ctx.dvelocities->cpu_data_p;
    auto &LJ_matrix = ctx.LJ_matrix->cpu_data_p;
    auto *excluded = ctx.excluded->cpu_data_p;
    auto *q_elscales = ctx.q_elscales->cpu_data_p;
    int ai, aj;
    real_t crg_i, crg_j;
    real_t elscale, scaling;
    bool bond23, bond14;
    coord_t da;
    real_t r2a, ra, r6a;
    real_t Vela, V_a, V_b;
    real_t dva;
    real_t ai_aii, aj_aii, ai_bii, aj_bii;
    real_t pair_a, pair_b, lookup_pair_a, lookup_pair_b, force_lj_term;

    for (int state = 0; state < ctx.n_lambdas; state++) {
        for (int qi = 0; qi < ctx.n_qatoms; qi++) {
            for (int qj = qi + 1; qj < ctx.n_qatoms; qj++) {
                ai = ctx.q_atoms[qi];
                aj = ctx.q_atoms[qj];

                crg_i = ctx.unified_ccharge(ai, state).charge;
                crg_j = ctx.unified_ccharge(aj, state).charge;

                bond23 = LJ_matrix[ai * ctx.n_atoms_solute + aj] == 3;
                bond14 = LJ_matrix[ai * ctx.n_atoms_solute + aj] == 1;

                if (bond23) continue;
                if (excluded[ai] || excluded[aj]) continue;

                scaling = bond14 ? ctx.topo.el14_scale : 1;

                elscale = 1;
                for (int k = 0; k < ctx.n_qelscales; k++) {
                    if (q_elscales[k + ctx.n_qelscales * state].qi == qi + 1 && q_elscales[k + ctx.n_qelscales * state].qj == qj + 1) {
                        elscale = q_elscales[k + ctx.n_qelscales * state].mu;
                    }
                }

                const catype_t& qi_type = ctx.unified_catype(ai, state);
                const catype_t& qj_type = ctx.unified_catype(aj, state);

                da.x = coords[aj].x - coords[ai].x;
                da.y = coords[aj].y - coords[ai].y;
                da.z = coords[aj].z - coords[ai].z;
                r2a = static_cast<real_t>(1.0) / (da.x * da.x + da.y * da.y + da.z * da.z);
                ra = static_cast<real_t>(std::sqrt(r2a));
                r6a = r2a * r2a * r2a;

                Vela = fortran_q_coulomb_energy(crg_i, crg_j, ra, ctx.topo.coulomb_constant, scaling * elscale);

                ai_aii = bond14 ? qi_type.aii_1_4 : qi_type.aii_normal;
                aj_aii = bond14 ? qj_type.aii_1_4 : qj_type.aii_normal;
                ai_bii = bond14 ? qi_type.bii_1_4 : qi_type.bii_normal;
                aj_bii = bond14 ? qj_type.bii_1_4 : qj_type.bii_normal;

                if (ctx.topo.vdw_rule == VDW_GEOMETRIC) {
                    calc_vdw_geometric(ai_aii, aj_aii, ai_bii, aj_bii, static_cast<real_t>(1.0), &pair_a, &pair_b);
                } else {
                    calc_vdw_arithmetic(ai_aii, aj_aii, ai_bii, aj_bii, static_cast<real_t>(1.0), &pair_a, &pair_b);
                }
                lookup_pair_a = pair_a;
                lookup_pair_b = pair_b;
                if (bond14 && ctx.q_softcore_use_max_potential) {
                    if (ctx.topo.vdw_rule == VDW_GEOMETRIC) {
                        calc_vdw_geometric(qi_type.aii_normal, qj_type.aii_normal, qi_type.bii_normal, qj_type.bii_normal, static_cast<real_t>(1.0), &lookup_pair_a, &lookup_pair_b);
                    } else {
                        calc_vdw_arithmetic(qi_type.aii_normal, qj_type.aii_normal, qi_type.bii_normal, qj_type.bii_normal, static_cast<real_t>(1.0), &lookup_pair_a, &lookup_pair_b);
                    }
                }
                const real_t softcore_alpha = q_softcore_pair_value(
                    ctx.q_softcore_value(qi, state),
                    ctx.q_softcore_value(qj, state),
                    ctx.q_softcore_use_max_potential);
                const real_t softcore_lookup = q_softcore_lookup_value(
                    softcore_alpha,
                    lookup_pair_a,
                    lookup_pair_b,
                    ctx.q_softcore_use_max_potential);
                q_softcore_lj(pair_a, pair_b, r6a, softcore_lookup, &V_a, &V_b, &force_lj_term);
                dva = r2a * (-Vela - force_lj_term) *
                      static_cast<real_t>(lambdas[state]);

                dvelocities[ai].x -= dva * da.x;
                dvelocities[ai].y -= dva * da.y;
                dvelocities[ai].z -= dva * da.z;

                dvelocities[aj].x += dva * da.x;
                dvelocities[aj].y += dva * da.y;
                dvelocities[aj].z += dva * da.z;

                ctx.EQ_nonbond_qq[state].Ucoul += static_cast<real_t>(Vela);
                ctx.EQ_nonbond_qq[state].Uvdw += static_cast<real_t>(V_a - V_b);
            }
        }
    }

#ifdef DEBUG
    printf("q-q: Ecoul = %f Evdw = %f\n", ctx.EQ_nonbond_qq[0].Ucoul, ctx.EQ_nonbond_qq[0].Uvdw);
#endif
}
