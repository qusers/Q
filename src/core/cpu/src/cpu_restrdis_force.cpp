#include "cpu_restrdis_force.h"

#include <math.h>

#include "context.h"

void calc_restrdis_forces() {
    auto& ctx = Context::instance();
    auto &coords = ctx.coords->cpu_data_p;
    auto &dvelocities = ctx.dvelocities->cpu_data_p;
    auto &restrdists = ctx.restrdists->cpu_data_p;
    auto *lambdas = ctx.lambdas->cpu_data_p;
    auto *EQ_restraint = ctx.EQ_restraint->cpu_data_p;

    int state, i, j;
    coord_t dr;
    real_t lambda, b, db, dv, ener;

    for (int ir = 0; ir < ctx.n_restrdists; ir++) {
        state = restrdists[ir].ipsi - 1;
        i = restrdists[ir].ai - 1;
        j = restrdists[ir].aj - 1;

        dr.x = coords[j].x - coords[i].x;
        dr.y = coords[j].y - coords[i].y;
        dr.z = coords[j].z - coords[i].z;

        if (restrdists[ir].ipsi != 0) {
            lambda = lambdas[state];
        } else {
            lambda = 1;
        }

        b = sqrt(pow(dr.x, 2) + pow(dr.y, 2) + pow(dr.z, 2));
        if (b < restrdists[ir].d1) {
            db = b - restrdists[ir].d1;
        } else if (b > restrdists[ir].d2) {
            db = b - restrdists[ir].d2;
        } else {
            continue;
        }

        ener = .5 * restrdists[ir].k * pow(db, 2);
        dv = lambda * restrdists[ir].k * db / b;

        dvelocities[j].x += dr.x * dv;
        dvelocities[j].y += dr.y * dv;
        dvelocities[j].z += dr.z * dv;
        dvelocities[i].x -= dr.x * dv;
        dvelocities[i].y -= dr.y * dv;
        dvelocities[i].z -= dr.z * dv;

        if (restrdists[ir].ipsi == 0) {
            for (int k = 0; k < ctx.n_lambdas; k++) {
                EQ_restraint[k].Urestr += ener;
            }
            if (ctx.n_lambdas == 0) {
                ctx.E_restraint.Upres += ener;
            }
        } else {
            EQ_restraint[state].Urestr += ener;
        }
    }
}
