#include "cpu_restrdis_force.h"

#include <math.h>

#include "context.h"

void calc_restrdis_forces() {
    auto& ctx = Context::instance();
    auto &coords = ctx.coords->cpu_data_p;
    auto &dvelocities = ctx.dvelocities->cpu_data_p;

    int state, i, j;
    coord_t dr;
    double lambda, b, db, dv, ener;

    for (int ir = 0; ir < ctx.n_restrdists; ir++) {
        state = ctx.restrdists[ir].ipsi - 1;
        i = ctx.restrdists[ir].ai - 1;
        j = ctx.restrdists[ir].aj - 1;

        dr.x = coords[j].x - coords[i].x;
        dr.y = coords[j].y - coords[i].y;
        dr.z = coords[j].z - coords[i].z;

        if (ctx.restrdists[ir].ipsi != 0) {
            lambda = ctx.lambdas[state];
        } else {
            lambda = 1;
        }

        b = sqrt(pow(dr.x, 2) + pow(dr.y, 2) + pow(dr.z, 2));
        if (b < ctx.restrdists[ir].d1) {
            db = b - ctx.restrdists[ir].d1;
        } else if (b > ctx.restrdists[ir].d2) {
            db = b - ctx.restrdists[ir].d2;
        } else {
            continue;
        }

        ener = .5 * ctx.restrdists[ir].k * pow(db, 2);
        dv = lambda * ctx.restrdists[ir].k * db / b;

        dvelocities[j].x += dr.x * dv;
        dvelocities[j].y += dr.y * dv;
        dvelocities[j].z += dr.z * dv;
        dvelocities[i].x -= dr.x * dv;
        dvelocities[i].y -= dr.y * dv;
        dvelocities[i].z -= dr.z * dv;

        if (ctx.restrdists[ir].ipsi == 0) {
            for (int k = 0; k < ctx.n_lambdas; k++) {
                ctx.EQ_restraint[k].Urestr += ener;
            }
            if (ctx.n_lambdas == 0) {
                ctx.E_restraint.Upres += ener;
            }
        } else {
            ctx.EQ_restraint[state].Urestr += ener;
        }
    }
}
