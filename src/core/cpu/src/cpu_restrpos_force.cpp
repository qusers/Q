#include "cpu_restrpos_force.h"

#include <math.h>

#include "context.h"

void calc_restrpos_forces() {
    auto& ctx = Context::instance();
    auto &coords = ctx.coords->cpu_data_p;
    auto &dvelocities = ctx.dvelocities->cpu_data_p;
    auto &restrspos = ctx.restrspos->cpu_data_p;
    auto *lambdas = ctx.lambdas->cpu_data_p;

    int state, i;
    coord_t dr;
    double lambda, ener, x2, y2, z2;

    for (int ir = 0; ir < ctx.n_restrspos; ir++) {
        state = restrspos[ir].ipsi - 1;
        i = restrspos[ir].a - 1;

        dr.x = coords[i].x - restrspos[ir].x.x;
        dr.y = coords[i].y - restrspos[ir].x.y;
        dr.z = coords[i].z - restrspos[ir].x.z;

        if (restrspos[ir].ipsi != 0) {
            lambda = lambdas[state];
        } else {
            lambda = 1;
        }

        x2 = pow(dr.x, 2);
        y2 = pow(dr.y, 2);
        z2 = pow(dr.z, 2);

        ener = .5 * restrspos[ir].k.x * x2 + .5 * restrspos[ir].k.y * y2 + .5 * restrspos[ir].k.z * z2;

        dvelocities[i].x += restrspos[ir].k.x * dr.x * lambda;
        dvelocities[i].y += restrspos[ir].k.y * dr.y * lambda;
        dvelocities[i].z += restrspos[ir].k.z * dr.z * lambda;

        if (restrspos[ir].ipsi == 0) {
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
