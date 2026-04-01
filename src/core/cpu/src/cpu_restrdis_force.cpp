#include "cpu_restrdis_force.h"

#include <math.h>

#include "context.h"

void calc_restrdis_forces() {
    auto& ctx = Context::instance();

    int state, i, j;
    coord_t dr;
    double lambda, b, db, dv, ener;

    for (int ir = 0; ir < ctx.n_restrdists; ir++) {
        state = ctx.restrdists[ir].ipsi - 1;
        i = ctx.restrdists[ir].ai - 1;
        j = ctx.restrdists[ir].aj - 1;

        dr.x = ctx.coords[j].x - ctx.coords[i].x;
        dr.y = ctx.coords[j].y - ctx.coords[i].y;
        dr.z = ctx.coords[j].z - ctx.coords[i].z;

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

        ctx.dvelocities[j].x += dr.x * dv;
        ctx.dvelocities[j].y += dr.y * dv;
        ctx.dvelocities[j].z += dr.z * dv;
        ctx.dvelocities[i].x -= dr.x * dv;
        ctx.dvelocities[i].y -= dr.y * dv;
        ctx.dvelocities[i].z -= dr.z * dv;

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
