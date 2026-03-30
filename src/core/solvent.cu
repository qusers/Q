#include "cpu/include/context.h"
#include "cpu/include/constants.h"

#include "system.h"
#include "solvent.h"
#include "vdw_rules.h"

#include <cmath>

namespace {

void init_water_parameters(Context &ctx) {
    const int oxygen_index = ctx.n_atoms_solute;
    const int hydrogen_index = oxygen_index + 1;

    const catype_t catype_ow = ctx.catypes[ctx.atypes[oxygen_index].code - 1];
    const ccharge_t ccharge_ow = ctx.ccharges[ctx.charges[oxygen_index].code - 1];
    const ccharge_t ccharge_hw = ctx.ccharges[ctx.charges[hydrogen_index].code - 1];

    if (ctx.topo.vdw_rule == VDW_GEOMETRIC) {
        ctx.A_OO = catype_ow.aii_normal * catype_ow.aii_normal;
        ctx.B_OO = catype_ow.bii_normal * catype_ow.bii_normal;
    } else {
        const double rstar_oo = 2.0 * catype_ow.aii_normal;
        const double eps_oo = catype_ow.bii_normal * catype_ow.bii_normal;
        const double r2_oo = rstar_oo * rstar_oo;
        const double r6_oo = r2_oo * r2_oo * r2_oo;
        ctx.A_OO = eps_oo * r6_oo * r6_oo;
        ctx.B_OO = 2.0 * eps_oo * r6_oo;
    }

    ctx.crg_ow = ccharge_ow.charge;
    ctx.crg_hw = ccharge_hw.charge;
}

inline double water_charge(const Context &ctx, int atom_offset) {
    return atom_offset == 0 ? ctx.crg_ow : ctx.crg_hw;
}

void accumulate_pair_force(Context &ctx,
                           int atom_i,
                           int atom_j,
                           double qi,
                           double qj,
                           bool include_vdw,
                           double vdw_a,
                           double vdw_b,
                           E_nonbonded_t &energy) {
    const double dx = ctx.coords[atom_j].x - ctx.coords[atom_i].x;
    const double dy = ctx.coords[atom_j].y - ctx.coords[atom_i].y;
    const double dz = ctx.coords[atom_j].z - ctx.coords[atom_i].z;

    const double r2inv = 1.0 / (dx * dx + dy * dy + dz * dz);
    const double rinv = std::sqrt(r2inv);
    const double ecoul = ctx.topo.coulomb_constant * qi * qj * rinv;

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

    ctx.dvelocities[atom_i].x -= scale * dx;
    ctx.dvelocities[atom_i].y -= scale * dy;
    ctx.dvelocities[atom_i].z -= scale * dz;

    ctx.dvelocities[atom_j].x += scale * dx;
    ctx.dvelocities[atom_j].y += scale * dy;
    ctx.dvelocities[atom_j].z += scale * dz;

    energy.Ucoul += ecoul;
    energy.Uvdw += evdw;
}

}  // namespace

void calc_nonbonded_ww_forces() {
    auto &ctx = Context::instance();
    if (ctx.n_waters <= 1) {
        return;
    }

    init_water_parameters(ctx);

    for (int water_i = 0; water_i < ctx.n_waters; ++water_i) {
        const int base_i = ctx.n_atoms_solute + 3 * water_i;
        for (int water_j = water_i + 1; water_j < ctx.n_waters; ++water_j) {
            const int base_j = ctx.n_atoms_solute + 3 * water_j;
            for (int atom_i = 0; atom_i < 3; ++atom_i) {
                for (int atom_j = 0; atom_j < 3; ++atom_j) {
                    accumulate_pair_force(ctx,
                                          base_i + atom_i,
                                          base_j + atom_j,
                                          water_charge(ctx, atom_i),
                                          water_charge(ctx, atom_j),
                                          atom_i == 0 && atom_j == 0,
                                          ctx.A_OO,
                                          ctx.B_OO,
                                          ctx.E_nonbond_ww);
                }
            }
        }
    }
}


void calc_nonbonded_pw_forces() {
    auto &ctx = Context::instance();
    if (ctx.n_waters == 0 || ctx.n_patoms == 0) {
        return;
    }

    for (int pi = 0; pi < ctx.n_patoms; ++pi) {
        const int atom_i = ctx.p_atoms[pi].a - 1;
        for (int atom_j = ctx.n_atoms_solute; atom_j < ctx.n_atoms; ++atom_j) {
            if (ctx.excluded[atom_i] || ctx.excluded[atom_j]) {
                continue;
            }

            const double qi = ctx.ccharges[ctx.charges[atom_i].code - 1].charge;
            const double qj = ctx.ccharges[ctx.charges[atom_j].code - 1].charge;

            const catype_t atom_i_type = ctx.catypes[ctx.atypes[atom_i].code - 1];
            const catype_t atom_j_type = ctx.catypes[ctx.atypes[atom_j].code - 1];

            double v_a = 0.0;
            double v_b = 0.0;
            const double dx = ctx.coords[atom_j].x - ctx.coords[atom_i].x;
            const double dy = ctx.coords[atom_j].y - ctx.coords[atom_i].y;
            const double dz = ctx.coords[atom_j].z - ctx.coords[atom_i].z;
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

            ctx.dvelocities[atom_i].x -= scale * dx;
            ctx.dvelocities[atom_i].y -= scale * dy;
            ctx.dvelocities[atom_i].z -= scale * dz;

            ctx.dvelocities[atom_j].x += scale * dx;
            ctx.dvelocities[atom_j].y += scale * dy;
            ctx.dvelocities[atom_j].z += scale * dz;

            ctx.E_nonbond_pw.Ucoul += ecoul;
            ctx.E_nonbond_pw.Uvdw += (v_a - v_b);
        }
    }
}


