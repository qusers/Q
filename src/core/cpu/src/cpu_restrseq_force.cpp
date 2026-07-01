#include "cpu_restrseq_force.h"

#include <math.h>

#include "context.h"
#include "cpu_force_accumulation.h"

void calc_restrseq_forces() {
    auto& ctx = Context::instance();
    auto& atypes = ctx.atypes->cpu_data_p;
    auto& catypes = ctx.catypes->cpu_data_p;
    auto& coords = ctx.coords->cpu_data_p;
    auto& dvelocities = ctx.dvelocities->cpu_data_p;
    auto& restrseqs = ctx.restrseqs->cpu_data_p;
    auto* heavy = ctx.heavy->cpu_data_p;

    energy_accum_t upres = 0;

    for (int s = 0; s < ctx.n_restrseqs; s++) {
        const double k = restrseqs[s].k;

        double dx = 0.0;
        double dy = 0.0;
        double dz = 0.0;
        int n_ctr = 0;
        double totmass = 0.0;

        if (restrseqs[s].to_center == 1) {
            for (int i = restrseqs[s].ai - 1; i < restrseqs[s].aj - 1; i++) {
                if (heavy[i] || restrseqs[s].ih) {
                    n_ctr++;
                    dx += coords[i].x - ctx.coords_init->cpu_data_p[i].x;
                    dy += coords[i].y - ctx.coords_init->cpu_data_p[i].y;
                    dz += coords[i].z - ctx.coords_init->cpu_data_p[i].z;
                }
            }

            if (n_ctr > 0) {
                const double inv_n_ctr = 1.0 / static_cast<double>(n_ctr);
                dx *= inv_n_ctr;
                dy *= inv_n_ctr;
                dz *= inv_n_ctr;
                const double r2 = dx * dx + dy * dy + dz * dz;
                const double ener = 0.5 * k * r2;
                add_energy(upres, ener);

                for (int i = restrseqs[s].ai - 1; i < restrseqs[s].aj - 1; i++) {
                    if (heavy[i] || restrseqs[s].ih) {
                        const double mass = catypes[atypes[i].code - 1].m;
                        const double tmp = 12.010;
                        add_force(dvelocities[i].x, k * dx * mass / tmp);
                        add_force(dvelocities[i].y, k * dy * mass / tmp);
                        add_force(dvelocities[i].z, k * dz * mass / tmp);
                    }
                }
            }
        } else if (restrseqs[s].to_center == 2) {
            for (int i = restrseqs[s].ai - 1; i < restrseqs[s].aj - 1; i++) {
                if (heavy[i] || restrseqs[s].ih) {
                    const double mass = catypes[atypes[i].code - 1].m;
                    totmass += mass;
                    dx += (coords[i].x - ctx.coords_init->cpu_data_p[i].x) * mass;
                    dy += (coords[i].y - ctx.coords_init->cpu_data_p[i].y) * mass;
                    dz += (coords[i].z - ctx.coords_init->cpu_data_p[i].z) * mass;
                }
            }

            if (totmass > 0) {
                dx /= totmass;
                dy /= totmass;
                dz /= totmass;
                const double r2 = dx * dx + dy * dy + dz * dz;
                const double ener = 0.5 * k * r2;
                add_energy(upres, ener);

                for (int i = restrseqs[s].ai - 1; i < restrseqs[s].aj - 1; i++) {
                    if (heavy[i] || restrseqs[s].ih) {
                        add_force(dvelocities[i].x, k * dx);
                        add_force(dvelocities[i].y, k * dy);
                        add_force(dvelocities[i].z, k * dz);
                    }
                }
            }
        } else {
            for (int i = restrseqs[s].ai - 1; i < restrseqs[s].aj - 1; i++) {
                if (heavy[i] || restrseqs[s].ih) {
                    const double dx_atom = coords[i].x - ctx.coords_init->cpu_data_p[i].x;
                    const double dy_atom = coords[i].y - ctx.coords_init->cpu_data_p[i].y;
                    const double dz_atom = coords[i].z - ctx.coords_init->cpu_data_p[i].z;

                    const double r2 = dx_atom * dx_atom + dy_atom * dy_atom + dz_atom * dz_atom;
                    const double ener = 0.5 * k * r2;
                    add_energy(upres, ener);

                    add_force(dvelocities[i].x, k * dx_atom);
                    add_force(dvelocities[i].y, k * dy_atom);
                    add_force(dvelocities[i].z, k * dz_atom);
                }
            }
        }
    }

    auto* energy = ctx.energy.host();
    energy[E_RESTR_PRES] += upres;
}
