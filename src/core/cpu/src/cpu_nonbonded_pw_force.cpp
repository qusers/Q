#include "cpu_nonbonded_pw_force.h"

#include <math.h>

#include "constants.h"
#include "context.h"
#include "vdw_rules.h"

void calc_nonbonded_pw_forces() {
    auto& ctx = Context::instance();
    auto &coords = ctx.coords->cpu_data_p;
    auto &dvelocities = ctx.dvelocities->cpu_data_p;
    auto *excluded = ctx.excluded->cpu_data_p;
    if (ctx.n_waters == 0 || ctx.n_patoms == 0) {
        return;
    }

    for (int pi = 0; pi < ctx.n_patoms; ++pi) {
        const int atom_i = ctx.p_atoms[pi];
        for (int atom_j = ctx.n_atoms_solute; atom_j < ctx.n_atoms; ++atom_j) {
            if (excluded[atom_i] || excluded[atom_j]) {
                continue;
            }

            const double qi = ctx.unified_ccharge(atom_i, 0).charge;
            const double qj = ctx.unified_ccharge(atom_j, 0).charge;

            const catype_t& atom_i_type = ctx.unified_catype(atom_i, 0);
            const catype_t& atom_j_type = ctx.unified_catype(atom_j, 0);

            double v_a = 0.0;
            double v_b = 0.0;
            const double dx = coords[atom_j].x - coords[atom_i].x;
            const double dy = coords[atom_j].y - coords[atom_i].y;
            const double dz = coords[atom_j].z - coords[atom_i].z;
            const double r2inv = 1.0 / (dx * dx + dy * dy + dz * dz);
            const double rinv = std::sqrt(r2inv);
            const double r6inv = r2inv * r2inv * r2inv;
            const double ecoul = ctx.topo.coulomb_constant * qi * qj * rinv;

            if (ctx.topo.vdw_rule == VDW_GEOMETRIC) {
                calc_vdw_geometric(atom_i_type.aii_normal,
                                   atom_j_type.aii_normal,
                                   atom_i_type.bii_normal,
                                   atom_j_type.bii_normal,
                                   r6inv,
                                   &v_a,
                                   &v_b);
            } else {
                calc_vdw_arithmetic(atom_i_type.aii_normal,
                                    atom_j_type.aii_normal,
                                    atom_i_type.bii_normal,
                                    atom_j_type.bii_normal,
                                    r6inv,
                                    &v_a,
                                    &v_b);
            }

            const double scale = r2inv * (-ecoul - 12.0 * v_a + 6.0 * v_b);

            dvelocities[atom_i].x -= scale * dx;
            dvelocities[atom_i].y -= scale * dy;
            dvelocities[atom_i].z -= scale * dz;

            dvelocities[atom_j].x += scale * dx;
            dvelocities[atom_j].y += scale * dy;
            dvelocities[atom_j].z += scale * dz;

            ctx.E_nonbond_pw.Ucoul += ecoul;
            ctx.E_nonbond_pw.Uvdw += (v_a - v_b);
        }
    }
}
