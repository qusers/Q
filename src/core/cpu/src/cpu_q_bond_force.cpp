#include "cpu_q_bond_force.h"

#include <math.h>

#include "context.h"

void calc_qbond_forces(int state) {
    auto& ctx = Context::instance();
    auto &coords = ctx.coords->cpu_data_p;
    auto &dvelocities = ctx.dvelocities->cpu_data_p;
    int ic;
    int ai, aj;
    double b, db, ener, dv;
    coord_t rij;

    for (int i = 0; i < ctx.n_qbonds; i++) {
        ic = ctx.q_bonds[i + ctx.n_qbonds * state].code;

        if (ic == 0) {
            continue;
        }

        ai = ctx.q_bonds[i + ctx.n_qbonds * state].ai - 1;
        aj = ctx.q_bonds[i + ctx.n_qbonds * state].aj - 1;

        rij.x = coords[aj].x - coords[ai].x;
        rij.y = coords[aj].y - coords[ai].y;
        rij.z = coords[aj].z - coords[ai].z;

        b = sqrt(pow(rij.x, 2) + pow(rij.y, 2) + pow(rij.z, 2));
        db = b - ctx.q_cbonds[ic].b0;

        ener = 0.5 * ctx.q_cbonds[ic].kb * pow(db, 2);
        ctx.EQ_bond[state].Ubond += ener;
        dv = db * ctx.q_cbonds[ic].kb * ctx.lambdas[state] / b;

        dvelocities[ai].x -= dv * rij.x;
        dvelocities[ai].y -= dv * rij.y;
        dvelocities[ai].z -= dv * rij.z;
        dvelocities[aj].x += dv * rij.x;
        dvelocities[aj].y += dv * rij.y;
        dvelocities[aj].z += dv * rij.z;
    }
}
