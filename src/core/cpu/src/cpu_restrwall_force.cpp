#include "cpu_restrwall_force.h"

#include <math.h>

#include "context.h"

void calc_restrwall_forces() {
    auto& ctx = Context::instance();

    double k, b, db, ener, dv, fexp;
    coord_t dr;

    for (int ir = 0; ir < ctx.n_restrwalls; ir++) {
        k = ctx.restrwalls[ir].k;
        for (int i = ctx.restrwalls[ir].ai - 1; i < ctx.restrwalls[ir].aj - 1; i++) {
            if (ctx.heavy[i] || ctx.restrwalls[ir].ih) {
                dr.x = ctx.coords[i].x - ctx.topo.solvent_center.x;
                dr.y = ctx.coords[i].y - ctx.topo.solvent_center.y;
                dr.z = ctx.coords[i].z - ctx.topo.solvent_center.z;

                b = sqrt(pow(dr.x, 2) + pow(dr.y, 2) + pow(dr.z, 2));
                db = b - ctx.restrwalls[ir].d;

                if (db > 0) {
                    ener = .5 * k * pow(db, 2) - ctx.restrwalls[ir].dMorse;
                    dv = k * db / b;
                } else {
                    fexp = exp(ctx.restrwalls[ir].aMorse * db);
                    ener = ctx.restrwalls[ir].dMorse * (fexp * fexp - 2 * fexp);
                    dv = -2 * ctx.restrwalls[ir].dMorse * ctx.restrwalls[ir].aMorse * (fexp - fexp * fexp) / b;
                }
                ctx.E_restraint.Upres += ener;

                ctx.dvelocities[i].x += dv * dr.x;
                ctx.dvelocities[i].y += dv * dr.y;
                ctx.dvelocities[i].z += dv * dr.z;
            }
        }
    }
}
