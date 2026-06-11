#include "cpu_nonbonded_pw_force.h"

#include <math.h>

#include "constants.h"
#include "context.h"
#include "vdw_rules.h"
#include "cpu_force_accumulation.h"

void calc_nonbonded_pw_forces() {
    auto& ctx = Context::instance();
    auto &coords = ctx.coords->cpu_data_p;
    auto &dvelocities = ctx.dvelocities->cpu_data_p;
    auto *excluded = ctx.excluded->cpu_data_p;
    energy_accum_t ucoul = 0;
    energy_accum_t uvdw = 0;
    if (ctx.n_waters == 0 || ctx.n_patoms == 0) {
        return;
    }

    for (int pi = 0; pi < ctx.n_patoms; ++pi) {
        const int atom_i = ctx.p_atoms[pi];
        for (int atom_j = ctx.n_atoms_solute; atom_j < ctx.n_atoms; ++atom_j) {
            if (excluded[atom_i] || excluded[atom_j]) {
                continue;
            }

            real_t qi = ctx.unified_ccharge(atom_i, 0).charge;
            real_t qj = ctx.unified_ccharge(atom_j, 0).charge;

            const catype_t& atom_i_type = ctx.unified_catype(atom_i, 0);
            const catype_t& atom_j_type = ctx.unified_catype(atom_j, 0);

            real_t v_a = 0.0;
            real_t v_b = 0.0;
            const real_t dx = coords[atom_j].x - coords[atom_i].x;
            const real_t dy = coords[atom_j].y - coords[atom_i].y;
            const real_t dz = coords[atom_j].z - coords[atom_i].z;
            const real_t r2inv = 1.0 / (dx * dx + dy * dy + dz * dz);
            const real_t rinv = (std::sqrt(r2inv));
            const real_t r6inv = r2inv * r2inv * r2inv;
            const real_t ecoul = qi * qj * rinv * ctx.topo.coulomb_constant;

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

            const real_t scale = r2inv * (-ecoul - static_cast<real_t>(12.0) * v_a + static_cast<real_t>(6.0) * v_b);

            add_force(dvelocities[atom_i].x, -scale * dx);
            add_force(dvelocities[atom_i].y, -scale * dy);
            add_force(dvelocities[atom_i].z, -scale * dz);

            add_force(dvelocities[atom_j].x, scale * dx);
            add_force(dvelocities[atom_j].y, scale * dy);
            add_force(dvelocities[atom_j].z, scale * dz);

            add_energy(ucoul, ecoul);
            add_energy(uvdw, v_a - v_b);
        }
    }

    ctx.E_nonbond_pw.Ucoul += energy_from_accum(ucoul);
    ctx.E_nonbond_pw.Uvdw += energy_from_accum(uvdw);
}
