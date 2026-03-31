#include "cpu_handler.h"

#include <time.h>

#include "context.h"
#include "cpu_bonded_force.h"
#include "cpu_leapfrog.h"
#include "cpu_nonbonded_pp_force.h"
#include "cpu_nonbonded_pw_force.h"
#include "cpu_polx_water_force.h"
#include "cpu_pshell_force.h"
#include "cpu_q_angle_force.h"
#include "cpu_q_bond_force.h"
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

namespace {

void calc_restraint_forces(int iteration, Context& ctx) {
    if (ctx.n_waters > 0) {
        calc_radix_w_forces();
        if (ctx.md.polarisation) {
            calc_polx_w_forces(iteration);
        }
    }

    calc_pshell_forces();
    calc_restrseq_forces();
    calc_restrdis_forces();
    calc_restrpos_forces();
    calc_restrang_forces();
    calc_restrwall_forces();
}

void calc_q_bonded_forces(Context& ctx) {
    for (int state = 0; state < ctx.n_lambdas; state++) {
        calc_qangle_forces(state);
        calc_qbond_forces(state);
        calc_qtorsion_forces(state);
    }
}

}  // namespace

void CpuHandler::initialize() {
}

void CpuHandler::shutdown() {
}

void CpuHandler::calc_internal_forces(int iteration) {

    calc_bonded_forces();
    calc_restraint_forces(iteration, ctx);
    calc_q_bonded_forces(ctx);
}

void CpuHandler::calc_nonbonded_forces() {
    calc_nonbonded_qp_forces();
    calc_nonbonded_pp_forces();
    calc_nonbonded_ww_forces();
    calc_nonbonded_pw_forces();
    calc_nonbonded_qw_forces();
    calc_nonbonded_qq_forces();
}

void CpuHandler::calc_temperature() {
    ::calc_temperature();
}

void CpuHandler::calc_leapfrog() {
    ::calc_leapfrog();
}
