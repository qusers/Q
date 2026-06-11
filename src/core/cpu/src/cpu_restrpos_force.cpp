#include "cpu_restrpos_force.h"

#include <math.h>
#include <vector>

#include "context.h"
#include "cpu_force_accumulation.h"

void calc_restrpos_forces() {
    auto& ctx = Context::instance();
    auto &coords = ctx.coords->cpu_data_p;
    auto &dvelocities = ctx.dvelocities->cpu_data_p;
    auto &restrspos = ctx.restrspos->cpu_data_p;
    auto *lambdas = ctx.lambdas->cpu_data_p;
    auto *EQ_restraint = ctx.EQ_restraint->cpu_data_p;

    int state, i;
    coord_t dr;
    real_t lambda, ener, x2, y2, z2;
    std::vector<energy_accum_t> urestr(ctx.n_lambdas, 0);
    energy_accum_t upres = 0;

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

        x2 = dr.x * dr.x;
        y2 = dr.y * dr.y;
        z2 = dr.z * dr.z;

        ener = .5 * restrspos[ir].k.x * x2 + .5 * restrspos[ir].k.y * y2 + .5 * restrspos[ir].k.z * z2;

        add_force(dvelocities[i].x, restrspos[ir].k.x * dr.x * lambda);
        add_force(dvelocities[i].y, restrspos[ir].k.y * dr.y * lambda);
        add_force(dvelocities[i].z, restrspos[ir].k.z * dr.z * lambda);

        if (restrspos[ir].ipsi == 0) {
            for (int k = 0; k < ctx.n_lambdas; k++) {
                add_energy(urestr[k], ener);
            }
            if (ctx.n_lambdas == 0) {
                add_energy(upres, ener);
            }
        } else {
            add_energy(urestr[state], ener);
        }
    }

    for (int state = 0; state < ctx.n_lambdas; state++) {
        EQ_restraint[state].Urestr += energy_from_accum(urestr[state]);
    }
    ctx.E_restraint.Upres += energy_from_accum(upres);
}
