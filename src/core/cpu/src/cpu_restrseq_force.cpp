#include "cpu_restrseq_force.h"

#include <math.h>

#include "context.h"
#include "cpu_force_accumulation.h"

void calc_restrseq_forces() {
    auto& ctx = Context::instance();
    auto &atypes = ctx.atypes->cpu_data_p;
    auto &catypes = ctx.catypes->cpu_data_p;
    auto &coords = ctx.coords->cpu_data_p;
    auto &dvelocities = ctx.dvelocities->cpu_data_p;
    auto &restrseqs = ctx.restrseqs->cpu_data_p;
    auto *heavy = ctx.heavy->cpu_data_p;

    real_t k, mass, totmass;
    coord_t dr;
    real_t r2, ener;
    energy_accum_t upres = 0;

    for (int s = 0; s < ctx.n_restrseqs; s++) {
        k = restrseqs[s].k;

        dr.x = 0;
        dr.y = 0;
        dr.z = 0;
        int n_ctr = 0;
        totmass = 0;

        if (restrseqs[s].to_center == 1) {
            for (int i = restrseqs[s].ai - 1; i < restrseqs[s].aj - 1; i++) {
                if (heavy[i] || restrseqs[s].ih) {
                    n_ctr++;
                    dr.x += (coords[i].x - ctx.coords_init->cpu_data_p[i].x);
                    dr.y += (coords[i].y - ctx.coords_init->cpu_data_p[i].y);
                    dr.z += (coords[i].z - ctx.coords_init->cpu_data_p[i].z);
                }
            }

            if (n_ctr > 0) {
                dr.x /= n_ctr;
                dr.y /= n_ctr;
                dr.z /= n_ctr;
                r2 = pow(dr.x, 2) + pow(dr.y, 2) + pow(dr.z, 2);
                ener = .5 * k * r2;
                add_energy(upres, ener);

                for (int i = restrseqs[s].ai - 1; i < restrseqs[s].aj - 1; i++) {
                    if (heavy[i] || restrseqs[s].ih) {
                        mass = catypes[atypes[i].code - 1].m;
                        const real_t tmp = 12.010;
                        add_force(dvelocities[i].x, k * dr.x * mass / tmp);
                        add_force(dvelocities[i].y, k * dr.y * mass / tmp);
                        add_force(dvelocities[i].z, k * dr.z * mass / tmp);
                    }
                }
            }
        } else if (restrseqs[s].to_center == 2) {
            for (int i = restrseqs[s].ai - 1; i < restrseqs[s].aj - 1; i++) {
                if (heavy[i] || restrseqs[s].ih) {
                    mass = catypes[atypes[i].code - 1].m;
                    totmass += mass;
                    dr.x += (coords[i].x - ctx.coords_init->cpu_data_p[i].x) * mass;
                    dr.y += (coords[i].y - ctx.coords_init->cpu_data_p[i].y) * mass;
                    dr.z += (coords[i].z - ctx.coords_init->cpu_data_p[i].z) * mass;
                }
            }

            if (totmass > 0) {
                dr.x /= totmass;
                dr.y /= totmass;
                dr.z /= totmass;
                r2 = pow(dr.x, 2) + pow(dr.y, 2) + pow(dr.z, 2);
                ener = .5 * k * r2;
                add_energy(upres, ener);

                for (int i = restrseqs[s].ai - 1; i < restrseqs[s].aj - 1; i++) {
                    if (heavy[i] || restrseqs[s].ih) {
                        add_force(dvelocities[i].x, k * dr.x);
                        add_force(dvelocities[i].y, k * dr.y);
                        add_force(dvelocities[i].z, k * dr.z);
                    }
                }
            }
        } else {
            for (int i = restrseqs[s].ai - 1; i < restrseqs[s].aj - 1; i++) {
                if (heavy[i] || restrseqs[s].ih) {
                    dr.x = coords[i].x - ctx.coords_init->cpu_data_p[i].x;
                    dr.y = coords[i].y - ctx.coords_init->cpu_data_p[i].y;
                    dr.z = coords[i].z - ctx.coords_init->cpu_data_p[i].z;

                    r2 = pow(dr.x, 2) + pow(dr.y, 2) + pow(dr.z, 2);
                    ener = .5 * k * r2;
                    add_energy(upres, ener);

                    add_force(dvelocities[i].x, k * dr.x);
                    add_force(dvelocities[i].y, k * dr.y);
                    add_force(dvelocities[i].z, k * dr.z);
                }
            }
        }
    }

    ctx.E_restraint.Upres += energy_from_accum(upres);
}
