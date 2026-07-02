#include "cpu_polx_water_force.h"
#include "cpu_force_accumulation.h"

#include <math.h>
#include <stdio.h>

#include "constants.h"
#include "context.h"

void calc_polx_w_forces(int iteration) {
    auto& ctx = Context::instance();
    auto &coords = ctx.coords->cpu_data_p;
    auto &dvelocities = ctx.dvelocities->cpu_data_p;
    auto &exclude = ctx.excluded->cpu_data_p;
    auto *wshells = ctx.wshells->cpu_data_p;

    int wi, imin, jw, ii, iis, jmin;
    double tmin;
    energy_accum_t upolx = 0;
    for (int is = 0; is < ctx.n_shells; is++) {
        wshells[is].n_inshell = 0;
    }

    for (int i = 0; i < ctx.n_waters; i++) {
        ctx.theta[i] = 0.0;
        ctx.theta0[i] = 0.0;

        wi = ctx.n_atoms_solute + 3 * i;
        if (exclude[wi]) continue;

        double rmu_x = coords[wi + 1].x +
                       coords[wi + 2].x -
                       2.0 * coords[wi].x;
        double rmu_y = coords[wi + 1].y +
                       coords[wi + 2].y -
                       2.0 * coords[wi].y;
        double rmu_z = coords[wi + 1].z +
                       coords[wi + 2].z -
                       2.0 * coords[wi].z;

        const double rm = sqrt(rmu_x * rmu_x + rmu_y * rmu_y + rmu_z * rmu_z);

        rmu_x /= rm;
        rmu_y /= rm;
        rmu_z /= rm;

        double rcu_x = coords[wi].x - ctx.topo.solvent_center.x;
        double rcu_y = coords[wi].y - ctx.topo.solvent_center.y;
        double rcu_z = coords[wi].z - ctx.topo.solvent_center.z;
        const double rc = sqrt(rcu_x * rcu_x + rcu_y * rcu_y + rcu_z * rcu_z);
        rcu_x /= rc;
        rcu_y /= rc;
        rcu_z /= rc;

        double cos_th = rmu_x * rcu_x + rmu_y * rcu_y + rmu_z * rcu_z;
        if (cos_th > 1.0) {
            cos_th = 1.0;
        }
        if (cos_th < -1.0) {
            cos_th = -1.0;
        }
        ctx.theta[i] = acos(cos_th);
        ctx.tdum[i] = ctx.theta[i];

        if (rc > wshells[ctx.n_shells - 1].router - wshells[ctx.n_shells - 1].dr) {
            for (iis = ctx.n_shells - 1; iis > 0; iis--) {
                if (rc <= wshells[iis].router) {
                    break;
                }
            }
            wshells[iis].n_inshell += 1;
            ctx.list_sh[wshells[iis].n_inshell - 1][iis] = i;
        }
    }
    const double PI = 4.0 * atan(1.0);
    for (int is = 0; is < ctx.n_shells; is++) {
        imin = 0;
        for (int il = 0; il < wshells[is].n_inshell; il++) {
            tmin = 2.0 * PI;
            for (int jl = 0; jl < wshells[is].n_inshell; jl++) {
                jw = ctx.list_sh[jl][is];
                if (ctx.tdum[jw] < tmin) {
                    jmin = jw;
                    tmin = ctx.theta[jw];
                }
            }
            ctx.nsort[imin][is] = jmin;
            imin++;
            ctx.tdum[jmin] = 99999.0;
        }
    }

    if (iteration != 0 && iteration % itdis_update == 0) {
        for (int is = 0; is < ctx.n_shells; is++) {
            printf("SHELL %d\n", is);
            wshells[is].avtheta /= itdis_update;
            wshells[is].avn_inshell /= itdis_update;
            wshells[is].theta_corr =
                wshells[is].theta_corr + wshells[is].avtheta - acos(wshells[is].cstb);
            printf("average theta = %f, average in shell = %f, theta_corr = %f\n",
                   wshells[is].avtheta * 180 / M_PI, wshells[is].avn_inshell,
                   wshells[is].theta_corr * 180 / M_PI);
            wshells[is].avtheta = 0.0;
            wshells[is].avn_inshell = 0.0;
        }
    }

    for (int is = 0; is < ctx.n_shells; is++) {
        if (wshells[is].n_inshell == 0) {
            continue;
        }

        double avtdum = 0.0;
        for (int il = 0; il < wshells[is].n_inshell; il++) {
            ii = ctx.nsort[il][is];
            const double arg = 1.0 + ((1.0 - 2.0 * (il + 1)) / wshells[is].n_inshell);
            ctx.theta0[il] = acos(arg);
            ctx.theta0[il] = ctx.theta0[il] - 3.0 * sin(ctx.theta0[il]) * wshells[is].cstb / 2.0;
            if (ctx.theta0[il] < 0.0) {
                ctx.theta0[il] = 0.0;
            }
            if (ctx.theta0[il] > PI) {
                ctx.theta0[il] = PI;
            }



            avtdum += ctx.theta[ii];
            const double dtheta = ctx.theta[ii] - ctx.theta0[il] + wshells[is].theta_corr;
            const double ener = 0.5 * ctx.md.polarisation_force * dtheta * dtheta;
            add_energy(upolx, ener);

            const double dv = ctx.md.polarisation_force * dtheta;

            wi = ctx.n_atoms_solute + 3 * ii;

            double rmu_x = coords[wi + 1].x +
                           coords[wi + 2].x -
                           2.0 * coords[wi].x;
            double rmu_y = coords[wi + 1].y +
                           coords[wi + 2].y -
                           2.0 * coords[wi].y;
            double rmu_z = coords[wi + 1].z +
                           coords[wi + 2].z -
                           2.0 * coords[wi].z;

            const double rm = sqrt(rmu_x * rmu_x + rmu_y * rmu_y + rmu_z * rmu_z);

            rmu_x /= rm;
            rmu_y /= rm;
            rmu_z /= rm;

            double rcu_x = coords[wi].x - ctx.topo.solvent_center.x;
            double rcu_y = coords[wi].y - ctx.topo.solvent_center.y;
            double rcu_z = coords[wi].z - ctx.topo.solvent_center.z;
            const double rc = sqrt(rcu_x * rcu_x + rcu_y * rcu_y + rcu_z * rcu_z);
            rcu_x /= rc;
            rcu_y /= rc;
            rcu_z /= rc;

            double cos_th = rmu_x * rcu_x + rmu_y * rcu_y + rmu_z * rcu_z;
            if (cos_th > 1.0) {
                cos_th = 1.0;
            }
            if (cos_th < -1.0) {
                cos_th = -1.0;
            }
            double f0 = sin(acos(cos_th));
            const double sin_epsilon = k_singular_sin_epsilon;
            if (fabs(f0) < sin_epsilon) {
                f0 = sin_epsilon;
            }
            f0 = -1.0 / f0;
            f0 *= dv;

            const double f1O_x = -2.0 * (rcu_x - rmu_x * cos_th) / rm;
            const double f1O_y = -2.0 * (rcu_y - rmu_y * cos_th) / rm;
            const double f1O_z = -2.0 * (rcu_z - rmu_z * cos_th) / rm;
            const double f1H1_x = (rcu_x - rmu_x * cos_th) / rm;
            const double f1H1_y = (rcu_y - rmu_y * cos_th) / rm;
            const double f1H1_z = (rcu_z - rmu_z * cos_th) / rm;
            const double f1H2_x = (rcu_x - rmu_x * cos_th) / rm;
            const double f1H2_y = (rcu_y - rmu_y * cos_th) / rm;
            const double f1H2_z = (rcu_z - rmu_z * cos_th) / rm;

            const double f2_x = (rmu_x - rcu_x * cos_th) / rc;
            const double f2_y = (rmu_y - rcu_y * cos_th) / rc;
            const double f2_z = (rmu_z - rcu_z * cos_th) / rc;

            add_force(dvelocities[wi].x, f0 * (f1O_x + f2_x));
            add_force(dvelocities[wi].y, f0 * (f1O_y + f2_y));
            add_force(dvelocities[wi].z, f0 * (f1O_z + f2_z));


            add_force(dvelocities[wi + 1].x, f0 * f1H1_x);
            add_force(dvelocities[wi + 1].y, f0 * f1H1_y);
            add_force(dvelocities[wi + 1].z, f0 * f1H1_z);

            add_force(dvelocities[wi + 2].x, f0 * f1H2_x);
            add_force(dvelocities[wi + 2].y, f0 * f1H2_y);
            add_force(dvelocities[wi + 2].z, f0 * f1H2_z);
        }

        wshells[is].avtheta += avtdum / wshells[is].n_inshell;
        wshells[is].avn_inshell += wshells[is].n_inshell;
    }

    auto* energy = ctx.energy.host();
    energy[E_RESTR_POLX] = upolx;
}
