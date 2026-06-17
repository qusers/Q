#include "cpu_nonbonded_ww_force.h"

#include <cmath>

#include "constants.h"
#include "context.h"
#include "vdw_rules.h"
#include "cpu_force_accumulation.h"

namespace {
void calc_oxygen_vdw_parameters(const Context& ctx, int oxygen_i, int oxygen_j, double* vdw_a, double* vdw_b) {
    const catype_t& oi_type = ctx.unified_catype(oxygen_i, 0);
    const catype_t& oj_type = ctx.unified_catype(oxygen_j, 0);
    if (ctx.topo.vdw_rule == VDW_GEOMETRIC) {
        *vdw_a = oi_type.aii_normal * oj_type.aii_normal;
        *vdw_b = oi_type.bii_normal * oj_type.bii_normal;
    } else {
        calc_vdw_arithmetic(
            oi_type.aii_normal, oj_type.aii_normal, oi_type.bii_normal, oj_type.bii_normal, 1.0, vdw_a, vdw_b);
    }
}
}  // namespace

void accumulate_pair_force(Context& ctx,
                           int atom_i,
                           int atom_j,
                           double qi,
                           double qj,
                           bool include_vdw,
                           double vdw_a,
                           double vdw_b,
                           energy_accum_t& ucoul,
                           energy_accum_t& uvdw) {
    auto &coords = ctx.coords->cpu_data_p;
    auto &dvelocities = ctx.dvelocities->cpu_data_p;
    const double dx = coords[atom_j].x - coords[atom_i].x;
    const double dy = coords[atom_j].y - coords[atom_i].y;
    const double dz = coords[atom_j].z - coords[atom_i].z;

    const double r2inv = 1.0 / (dx * dx + dy * dy + dz * dz);
    const double rinv = sqrt(r2inv);
    const double ecoul = qi * qj * rinv * ctx.topo.coulomb_constant;

    double evdw = 0.0;
    double dva = -ecoul;
    if (include_vdw) {
        const double r6inv = r2inv * r2inv * r2inv;
        const double v_a = vdw_a * r6inv * r6inv;
        const double v_b = vdw_b * r6inv;
        evdw = v_a - v_b;
        dva -= 12.0 * v_a - 6.0 * v_b;
    }

    const double scale = r2inv * dva;

    add_force(dvelocities[atom_i].x, -scale * dx);
    add_force(dvelocities[atom_i].y, -scale * dy);
    add_force(dvelocities[atom_i].z, -scale * dz);

    add_force(dvelocities[atom_j].x, scale * dx);
    add_force(dvelocities[atom_j].y, scale * dy);
    add_force(dvelocities[atom_j].z, scale * dz);

    add_energy(ucoul, ecoul);
    add_energy(uvdw, evdw);
}

void calc_nonbonded_ww_forces() {
    auto& ctx = Context::instance();
    if (ctx.n_waters <= 1) {
        return;
    }
    energy_accum_t ucoul = 0;
    energy_accum_t uvdw = 0;

    for (int water_i = 0; water_i < ctx.n_waters; ++water_i) {
        const int base_i = ctx.n_atoms_solute + 3 * water_i;
        for (int water_j = water_i + 1; water_j < ctx.n_waters; ++water_j) {
            const int base_j = ctx.n_atoms_solute + 3 * water_j;
            double oxygen_vdw_a = 0.0;
            double oxygen_vdw_b = 0.0;
            calc_oxygen_vdw_parameters(ctx, base_i, base_j, &oxygen_vdw_a, &oxygen_vdw_b);
            for (int atom_i = 0; atom_i < 3; ++atom_i) {
                for (int atom_j = 0; atom_j < 3; ++atom_j) {
                    accumulate_pair_force(ctx,
                                          base_i + atom_i,
                                          base_j + atom_j,
                                          ctx.unified_ccharge(base_i + atom_i, 0).charge,
                                          ctx.unified_ccharge(base_j + atom_j, 0).charge,
                                          atom_i == 0 && atom_j == 0,
                                          oxygen_vdw_a,
                                          oxygen_vdw_b,
                                          ucoul,
                                          uvdw);
                }
            }
        }
    }
    ctx.E_nonbond_ww.Ucoul += energy_from_accum(ucoul);
    ctx.E_nonbond_ww.Uvdw += energy_from_accum(uvdw);
}
