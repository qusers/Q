#include "cpu_nonbonded_ww_force.h"

#include "constants.h"
#include "context.h"
#include "vdw_rules.h"

void init_water_parameters(Context& ctx) {
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

inline double water_charge(const Context& ctx, int atom_offset) {
    return atom_offset == 0 ? ctx.crg_ow : ctx.crg_hw;
}

void accumulate_pair_force(Context& ctx,
                           int atom_i,
                           int atom_j,
                           double qi,
                           double qj,
                           bool include_vdw,
                           double vdw_a,
                           double vdw_b,
                           E_nonbonded_t& energy) {
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

void calc_nonbonded_ww_forces() {
    auto& ctx = Context::instance();
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