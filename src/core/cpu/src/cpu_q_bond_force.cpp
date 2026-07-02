#include "cpu_q_bond_force.h"
#include "cpu_force_accumulation.h"

#include <math.h>

#include "context.h"

void calc_qbond_forces(int state) {
    auto& ctx = Context::instance();
    auto &coords = ctx.coords->cpu_data_p;
    auto &dvelocities = ctx.dvelocities->cpu_data_p;
    auto *lambdas = ctx.lambdas->cpu_data_p;
    int ic;
    int ai, aj;
    energy_accum_t bond = 0;

    for (int i = 0; i < ctx.n_qbonds; i++) {
        ic = ctx.q_bonds[i + ctx.n_qbonds * state].code;

        if (ic == 0) {
            continue;
        }

        ai = ctx.q_bonds[i + ctx.n_qbonds * state].ai - 1;
        aj = ctx.q_bonds[i + ctx.n_qbonds * state].aj - 1;

        const double dx = coords[aj].x - coords[ai].x;
        const double dy = coords[aj].y - coords[ai].y;
        const double dz = coords[aj].z - coords[ai].z;

        const double b = sqrt(dx * dx + dy * dy + dz * dz);
        const double db = b - ctx.q_cbonds[ic].b0;

        const double ener = 0.5 * ctx.q_cbonds[ic].kb * db * db;
        add_energy(bond, ener);
        const double dv = db * ctx.q_cbonds[ic].kb * lambdas[state] / b;

        add_force(dvelocities[ai].x, -dv * dx);
        add_force(dvelocities[ai].y, -dv * dy);
        add_force(dvelocities[ai].z, -dv * dz);
        add_force(dvelocities[aj].x, dv * dx);
        add_force(dvelocities[aj].y, dv * dy);
        add_force(dvelocities[aj].z, dv * dz);
    }

    ctx.energy.host()[ctx.energy.eq_index(ENERGY_FIXED_COUNT, state, EQ_BOND_BOND)] += bond;
}
