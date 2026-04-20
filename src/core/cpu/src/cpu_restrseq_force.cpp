#include "cpu_restrseq_force.h"

#include <math.h>

#include "context.h"

void calc_restrseq_forces() {
    auto& ctx = Context::instance();

    double k, mass, totmass;
    coord_t dr;
    double r2, ener;

    for (int s = 0; s < ctx.n_restrseqs; s++) {
        k = ctx.restrseqs[s].k;

        dr.x = 0;
        dr.y = 0;
        dr.z = 0;
        int n_ctr = 0;
        totmass = 0;

        if (ctx.restrseqs[s].to_center == 1) {
            for (int i = ctx.restrseqs[s].ai - 1; i < ctx.restrseqs[s].aj - 1; i++) {
                if (ctx.heavy[i] || ctx.restrseqs[s].ih) {
                    n_ctr++;
                    dr.x += (ctx.coords[i].x - ctx.coords_init->cpu_data_p[i].x);
                    dr.y += (ctx.coords[i].y - ctx.coords_init->cpu_data_p[i].y);
                    dr.z += (ctx.coords[i].z - ctx.coords_init->cpu_data_p[i].z);
                }
            }

            if (n_ctr > 0) {
                dr.x /= n_ctr;
                dr.y /= n_ctr;
                dr.z /= n_ctr;
                r2 = pow(dr.x, 2) + pow(dr.y, 2) + pow(dr.z, 2);
                ener = .5 * k * r2;
                ctx.E_restraint.Upres += ener;

                for (int i = ctx.restrseqs[s].ai - 1; i < ctx.restrseqs[s].aj - 1; i++) {
                    if (ctx.heavy[i] || ctx.restrseqs[s].ih) {
                        mass = ctx.catypes[ctx.atypes[i].code - 1].m;
                        ctx.dvelocities[i].x += (k * dr.x * mass / 12.010);
                        ctx.dvelocities[i].y += (k * dr.y * mass / 12.010);
                        ctx.dvelocities[i].z += (k * dr.z * mass / 12.010);
                    }
                }
            }
        } else if (ctx.restrseqs[s].to_center == 2) {
            for (int i = ctx.restrseqs[s].ai - 1; i < ctx.restrseqs[s].aj - 1; i++) {
                if (ctx.heavy[i] || ctx.restrseqs[s].ih) {
                    mass = ctx.catypes[ctx.atypes[i].code - 1].m;
                    totmass += mass;
                    dr.x += (ctx.coords[i].x - ctx.coords_init->cpu_data_p[i].x) * mass;
                    dr.y += (ctx.coords[i].y - ctx.coords_init->cpu_data_p[i].y) * mass;
                    dr.z += (ctx.coords[i].z - ctx.coords_init->cpu_data_p[i].z) * mass;
                }
            }

            if (totmass > 0) {
                dr.x /= totmass;
                dr.y /= totmass;
                dr.z /= totmass;
                r2 = pow(dr.x, 2) + pow(dr.y, 2) + pow(dr.z, 2);
                ener = .5 * k * r2;
                ctx.E_restraint.Upres += ener;

                for (int i = ctx.restrseqs[s].ai - 1; i < ctx.restrseqs[s].aj - 1; i++) {
                    if (ctx.heavy[i] || ctx.restrseqs[s].ih) {
                        ctx.dvelocities[i].x += k * dr.x;
                        ctx.dvelocities[i].y += k * dr.y;
                        ctx.dvelocities[i].z += k * dr.z;
                    }
                }
            }
        } else {
            for (int i = ctx.restrseqs[s].ai - 1; i < ctx.restrseqs[s].aj - 1; i++) {
                if (ctx.heavy[i] || ctx.restrseqs[s].ih) {
                    dr.x = ctx.coords[i].x - ctx.coords_init->cpu_data_p[i].x;
                    dr.y = ctx.coords[i].y - ctx.coords_init->cpu_data_p[i].y;
                    dr.z = ctx.coords[i].z - ctx.coords_init->cpu_data_p[i].z;

                    r2 = pow(dr.x, 2) + pow(dr.y, 2) + pow(dr.z, 2);
                    ener = .5 * k * r2;
                    ctx.E_restraint.Upres += ener;

                    ctx.dvelocities[i].x += k * dr.x;
                    ctx.dvelocities[i].y += k * dr.y;
                    ctx.dvelocities[i].z += k * dr.z;
                }
            }
        }
    }
}
