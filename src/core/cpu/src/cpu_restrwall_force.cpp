#include "cpu_restrwall_force.h"

#include <math.h>

#include "context.h"
#include "cpu_force_accumulation.h"

void calc_restrwall_forces() {
    auto& ctx = Context::instance();
    auto& coords = ctx.coords->cpu_data_p;
    auto& dvelocities = ctx.dvelocities->cpu_data_p;
    auto& restrwalls = ctx.restrwalls->cpu_data_p;
    auto* heavy = ctx.heavy->cpu_data_p;

    energy_accum_t upres = 0;

    for (int ir = 0; ir < ctx.n_restrwalls; ir++) {
        const double k = restrwalls[ir].k;
        for (int i = restrwalls[ir].ai - 1; i < restrwalls[ir].aj - 1; i++) {
            if (heavy[i] || restrwalls[ir].ih) {
                const double dx = coords[i].x - ctx.topo.solvent_center.x;
                const double dy = coords[i].y - ctx.topo.solvent_center.y;
                const double dz = coords[i].z - ctx.topo.solvent_center.z;

                const double b = sqrt(dx * dx + dy * dy + dz * dz);
                const double db = b - restrwalls[ir].d;

                double ener;
                double dv;
                if (db > 0.0) {
                    ener = 0.5 * k * db * db - restrwalls[ir].dMorse;
                    dv = k * db / b;
                } else {
                    const double fexp = exp(restrwalls[ir].aMorse * db);
                    ener = restrwalls[ir].dMorse * (fexp * fexp - 2.0 * fexp);
                    dv = -2.0 * restrwalls[ir].dMorse * restrwalls[ir].aMorse * (fexp - fexp * fexp) / b;
                }
                add_energy(upres, ener);

                add_force(dvelocities[i].x, dv * dx);
                add_force(dvelocities[i].y, dv * dy);
                add_force(dvelocities[i].z, dv * dz);
            }
        }
    }

    auto* energy = ctx.energy.host();
    energy[EQ_RESTR_URESTR] += upres;
}
