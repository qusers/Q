#include "cpu_bonded_force.h"

#include "context.h"
#include "cpu_angle_force.h"
#include "cpu_bond_force.h"
#include "cpu_improper2_force.h"
#include "cpu_torsion_force.h"

void calc_bonded_forces() {
    auto& ctx = Context::instance();
    ctx.E_bond_p.Ubond = calc_bond_forces(0, ctx.n_bonds_solute);
    ctx.E_bond_w.Ubond = calc_bond_forces(ctx.n_bonds_solute, ctx.n_bonds);


    ctx.E_bond_p.Uangle = calc_angle_forces(0, ctx.n_angles_solute);
    ctx.E_bond_w.Uangle = calc_angle_forces(ctx.n_angles_solute, ctx.n_angles);

    ctx.E_bond_p.Utor = calc_torsion_forces(0, ctx.n_torsions_solute);
    ctx.E_bond_w.Utor = calc_torsion_forces(ctx.n_torsions_solute, ctx.n_torsions);

    ctx.E_bond_p.Uimp = calc_improper2_forces(0, ctx.n_impropers_solute);
    ctx.E_bond_w.Uimp = calc_improper2_forces(ctx.n_impropers_solute, ctx.n_impropers);
}
