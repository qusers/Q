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
    std::vector<energy_accum_t> urestr(ctx.n_lambdas, 0);
    energy_accum_t upres = 0;

    for (int ir = 0; ir < ctx.n_restrspos; ir++) {
        state = restrspos[ir].ipsi - 1;
        i = restrspos[ir].a - 1;

        const double dx = static_cast<double>(coords[i].x) - static_cast<double>(restrspos[ir].x.x);
        const double dy = static_cast<double>(coords[i].y) - static_cast<double>(restrspos[ir].x.y);
        const double dz = static_cast<double>(coords[i].z) - static_cast<double>(restrspos[ir].x.z);

        double lambda;
        if (restrspos[ir].ipsi != 0) {
            lambda = static_cast<double>(lambdas[state]);
        } else {
            lambda = 1.0;
        }

        const double kx = static_cast<double>(restrspos[ir].k.x);
        const double ky = static_cast<double>(restrspos[ir].k.y);
        const double kz = static_cast<double>(restrspos[ir].k.z);

        const double ener = 0.5 * kx * dx * dx + 0.5 * ky * dy * dy + 0.5 * kz * dz * dz;

        add_force(dvelocities[i].x, kx * dx * lambda);
        add_force(dvelocities[i].y, ky * dy * lambda);
        add_force(dvelocities[i].z, kz * dz * lambda);

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
