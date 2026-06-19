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
    std::vector<energy_accum_t> urestr(ctx.n_lambdas, 0);
    energy_accum_t upres = 0;

    for (int ir = 0; ir < ctx.n_restrdists; ir++) {
        state = restrdists[ir].ipsi - 1;
        i = restrdists[ir].ai - 1;
        j = restrdists[ir].aj - 1;

        const double dx = coords[j].x - coords[i].x;
        const double dy = coords[j].y - coords[i].y;
        const double dz = coords[j].z - coords[i].z;

        double lambda;
        if (restrdists[ir].ipsi != 0) {
            lambda = lambdas[state];
        } else {
            lambda = 1.0;
        }

        const double b = sqrt(dx * dx + dy * dy + dz * dz);
        double db;
        if (b < restrdists[ir].d1) {
            db = b - restrdists[ir].d1;
        } else if (b > restrdists[ir].d2) {
            db = b - restrdists[ir].d2;
        } else {
            continue;
        }

        const double ener = 0.5 * restrdists[ir].k * db * db;
        const double dv = lambda * restrdists[ir].k * db / b;

        add_force(dvelocities[j].x, dx * dv);
        add_force(dvelocities[j].y, dy * dv);
        add_force(dvelocities[j].z, dz * dv);
        add_force(dvelocities[i].x, -dx * dv);
        add_force(dvelocities[i].y, -dy * dv);
        add_force(dvelocities[i].z, -dz * dv);

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
