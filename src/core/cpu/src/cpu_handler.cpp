#include "cpu_handler.h"

#include <time.h>

#include "context.h"
#include "constants.h"
#include "cpu_angle_force.h"
#include "cpu_bond_force.h"
#include "cpu_bonded_force.h"
#include "cpu_improper2_force.h"
#include "cpu_leapfrog.h"
#include "cpu_nonbonded_pp_force.h"
#include "cpu_nonbonded_pw_force.h"
#include "cpu_polx_water_force.h"
#include "cpu_pshell_force.h"
#include "cpu_q_angle_force.h"
#include "cpu_q_bond_force.h"
#include "cpu_q_improper_force.h"
#include "cpu_nonbonded_qp_force.h"
#include "cpu_nonbonded_qq_force.h"
#include "cpu_nonbonded_qw_force.h"
#include "cpu_nonbonded_ww_force.h"
#include "cpu_q_torsion_force.h"
#include "cpu_radix_water_force.h"
#include "cpu_restrang_force.h"
#include "cpu_restrdis_force.h"
#include "cpu_restrpos_force.h"
#include "cpu_restrseq_force.h"
#include "cpu_restrwall_force.h"
#include "cpu_temperature.h"
#include "cpu_torsion_force.h"

namespace {

void calc_restraint_forces(int iteration, Context& ctx) {
    calc_pshell_forces();
    calc_restrseq_forces();
    calc_restrpos_forces();
    calc_restrdis_forces();
    calc_restrang_forces();
    calc_restrwall_forces();

    if (ctx.n_waters > 0) {
        calc_radix_w_forces();
        if (ctx.md.polarisation) {
            calc_polx_w_forces(iteration);
        }
    }
}

void calc_q_bonded_forces(Context& ctx) {
    for (int state = 0; state < ctx.n_lambdas; state++) {
        calc_qbond_forces(state);
        calc_qangle_forces(state);
        calc_qtorsion_forces(state);
        calc_qimproper_forces(state);
    }
}

}  // namespace

void CpuHandler::initialize_backend() {
}

void CpuHandler::shutdown() {
}

void CpuHandler::calc_internal_forces(int iteration) {
    calc_nonbonded_pp_forces();
    dump_force_group("01a_pp", iteration);
    if (ctx.topo.vdw_rule == VDW_ARITHMETIC) {
        calc_nonbonded_qp_forces();
        dump_force_group("01c_qp", iteration);
        calc_nonbonded_pw_forces();
        dump_force_group("01b_pw", iteration);
        calc_nonbonded_qw_forces();
        dump_force_group("01e_qw", iteration);
        calc_nonbonded_ww_forces();
        dump_force_group("01d_ww", iteration);
    } else {
        calc_nonbonded_pw_forces();
        dump_force_group("01b_pw", iteration);
        calc_nonbonded_qp_forces();
        dump_force_group("01c_qp", iteration);
        calc_nonbonded_ww_forces();
        dump_force_group("01d_ww", iteration);
        calc_nonbonded_qw_forces();
        dump_force_group("01e_qw", iteration);
    }
    dump_force_group("01_classical_nonbonds", iteration);

    ctx.E_bond_p.Ubond = calc_bond_forces(0, ctx.n_bonds_solute);
    ctx.E_bond_w.Ubond = calc_bond_forces(ctx.n_bonds_solute, ctx.n_bonds);
    dump_force_group("02a_bonds", iteration);

    ctx.E_bond_p.Uangle = calc_angle_forces(0, ctx.n_angles_solute);
    ctx.E_bond_w.Uangle = calc_angle_forces(ctx.n_angles_solute, ctx.n_angles);
    dump_force_group("02b_angles", iteration);

    ctx.E_bond_p.Utor = calc_torsion_forces(0, ctx.n_torsions_solute);
    ctx.E_bond_w.Utor = calc_torsion_forces(ctx.n_torsions_solute, ctx.n_torsions);
    dump_force_group("02c_torsions", iteration);

    ctx.E_bond_p.Uimp = calc_improper2_forces(0, ctx.n_impropers_solute);
    ctx.E_bond_w.Uimp = calc_improper2_forces(ctx.n_impropers_solute, ctx.n_impropers);
    dump_force_group("02d_impropers", iteration);
    dump_force_group("02_classical_bonded", iteration);
    calc_restraint_forces(iteration, ctx);
    dump_force_group("03_restraints", iteration);
    calc_nonbonded_qq_forces();
    dump_force_group("04_q_q_nonbonds", iteration);
    calc_q_bonded_forces(ctx);
    dump_force_group("05_q_bonded", iteration);
}

void CpuHandler::calc_nonbonded_forces() {
}

void CpuHandler::calc_temperature() {
    ::calc_temperature();
}

void CpuHandler::calc_leapfrog() {
    ::calc_leapfrog();
}
