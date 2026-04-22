#include "cpu_restrwall_force.h"

#include <math.h>

#include "context.h"

void calc_restrwall_forces() {
    auto& ctx = Context::instance();
    auto &coords = ctx.coords->cpu_data_p;
    auto &dvelocities = ctx.dvelocities->cpu_data_p;
    auto &restrwalls = ctx.restrwalls->cpu_data_p;
    auto *heavy = ctx.heavy->cpu_data_p;

    double k, b, db, ener, dv, fexp;
    coord_t dr;

    for (int ir = 0; ir < ctx.n_restrwalls; ir++) {
        k = restrwalls[ir].k;
        for (int i = restrwalls[ir].ai - 1; i < restrwalls[ir].aj - 1; i++) {
            if (heavy[i] || restrwalls[ir].ih) {
                dr.x = coords[i].x - ctx.topo.solvent_center.x;
                dr.y = coords[i].y - ctx.topo.solvent_center.y;
                dr.z = coords[i].z - ctx.topo.solvent_center.z;

                b = sqrt(pow(dr.x, 2) + pow(dr.y, 2) + pow(dr.z, 2));
                db = b - restrwalls[ir].d;

                if (db > 0) {
                    ener = .5 * k * pow(db, 2) - restrwalls[ir].dMorse;
                    dv = k * db / b;
                } else {
                    fexp = exp(restrwalls[ir].aMorse * db);
                    ener = restrwalls[ir].dMorse * (fexp * fexp - 2 * fexp);
                    dv = -2 * restrwalls[ir].dMorse * restrwalls[ir].aMorse * (fexp - fexp * fexp) / b;
                }
                ctx.E_restraint.Upres += ener;

                dvelocities[i].x += dv * dr.x;
                dvelocities[i].y += dv * dr.y;
                dvelocities[i].z += dv * dr.z;
            }
        }
    }
}
