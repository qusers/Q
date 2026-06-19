#include "cpu_nonbonded_pp_force.h"

#include <math.h>

#include "constants.h"
#include "context.h"
#include "vdw_rules.h"
#include "cpu_force_accumulation.h"

void calc_nonbonded_pp_forces() {
    auto& ctx = Context::instance();
    auto& coords = ctx.coords->cpu_data_p;
    auto& dvelocities = ctx.dvelocities->cpu_data_p;
    auto& LJ_matrix = ctx.LJ_matrix->cpu_data_p;
    auto* excluded = ctx.excluded->cpu_data_p;
    bool bond14, bond23;
    double scaling;
    coord_t da;
    double r2a, ra, r6a;
    double V_a, V_b;
    double crg_i, crg_j;
    double ai_aii, aj_aii, ai_bii, aj_bii;
    int i, j;
    energy_accum_t ucoul = 0;
    energy_accum_t uvdw = 0;

    for (int pi = 0; pi < ctx.n_patoms; pi++) {
        for (int pj = pi + 1; pj < ctx.n_patoms; pj++) {
            i = ctx.p_atoms[pi];
            j = ctx.p_atoms[pj];
            bond23 = LJ_matrix[i * ctx.n_atoms_solute + j] == 3;
            bond14 = LJ_matrix[i * ctx.n_atoms_solute + j] == 1;

            if (bond23) continue;
            if (excluded[i] || excluded[j]) continue;

            scaling = bond14 ? ctx.topo.el14_scale : 1;

            crg_i = ctx.unified_ccharge(i, 0).charge;
            crg_j = ctx.unified_ccharge(j, 0).charge;

            const catype_t& ai_type = ctx.unified_catype(i, 0);
            const catype_t& aj_type = ctx.unified_catype(j, 0);

            da.x = coords[j].x - coords[i].x;
            da.y = coords[j].y - coords[i].y;
            da.z = coords[j].z - coords[i].z;
            r2a = 1.0 / (da.x * da.x + da.y * da.y + da.z * da.z);
            ra = sqrt(r2a);
            r6a = r2a * r2a * r2a;

            const double Vela = crg_i * crg_j * ra * scaling * ctx.topo.coulomb_constant;

            ai_aii = bond14 ? ai_type.aii_1_4 : ai_type.aii_normal;
            aj_aii = bond14 ? aj_type.aii_1_4 : aj_type.aii_normal;
            ai_bii = bond14 ? ai_type.bii_1_4 : ai_type.bii_normal;
            aj_bii = bond14 ? aj_type.bii_1_4 : aj_type.bii_normal;

            if (ctx.topo.vdw_rule == VDW_GEOMETRIC) {
                calc_vdw_geometric(ai_aii, aj_aii, ai_bii, aj_bii, r6a, &V_a, &V_b);
            } else {
                calc_vdw_arithmetic(ai_aii, aj_aii, ai_bii, aj_bii, r6a, &V_a, &V_b);
            }
            const double dva = r2a * (-Vela - 12.0 * V_a + 6.0 * V_b);

            add_force(dvelocities[i].x, -dva * da.x);
            add_force(dvelocities[i].y, -dva * da.y);
            add_force(dvelocities[i].z, -dva * da.z);

            add_force(dvelocities[j].x, dva * da.x);
            add_force(dvelocities[j].y, dva * da.y);
            add_force(dvelocities[j].z, dva * da.z);

            add_energy(ucoul, Vela);
            add_energy(uvdw, V_a - V_b);
        }
    }

    ctx.E_nonbond_pp.Ucoul += energy_from_accum(ucoul);
    ctx.E_nonbond_pp.Uvdw += energy_from_accum(uvdw);
}
