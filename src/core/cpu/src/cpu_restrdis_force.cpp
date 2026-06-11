#include "cpu_restrdis_force.h"

#include <math.h>
#include <vector>

#include "context.h"
#include "cpu_force_accumulation.h"

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
    std::vector<energy_accum_t> urestr(ctx.n_lambdas, 0);
    energy_accum_t upres = 0;

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

        add_force(dvelocities[j].x, dr.x * dv);
        add_force(dvelocities[j].y, dr.y * dv);
        add_force(dvelocities[j].z, dr.z * dv);
        add_force(dvelocities[i].x, -dr.x * dv);
        add_force(dvelocities[i].y, -dr.y * dv);
        add_force(dvelocities[i].z, -dr.z * dv);

        if (restrdists[ir].ipsi == 0) {
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
