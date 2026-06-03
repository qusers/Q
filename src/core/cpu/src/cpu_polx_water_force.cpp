#include "cpu_polx_water_force.h"
#include "debug.h"

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
    real_t tmin;
    coord_t rmu, rcu, f1O, f1H1, f1H2, f2;
    real_t rm, rc;
    real_t cos_th;
    real_t avtdum, arg, f0, dv;
    real_t ener;
    std::ostringstream ss;
    ss << std::fixed << std::setprecision(8);
    for (int is = 0; is < ctx.n_shells; is++) {
        wshells[is].n_inshell = 0;
    }

    for (int i = 0; i < ctx.n_waters; i++) {
        ctx.theta[i] = 0;
        ctx.theta0[i] = 0;

        wi = ctx.n_atoms_solute + 3 * i;
        if (exclude[wi]) continue;

        rmu.x = coords[wi + 1].x + coords[wi + 2].x - 2 * coords[wi].x;
        rmu.y = coords[wi + 1].y + coords[wi + 2].y - 2 * coords[wi].y;
        rmu.z = coords[wi + 1].z + coords[wi + 2].z - 2 * coords[wi].z;

        rm = sqrt(rmu.x * rmu.x + rmu.y * rmu.y + rmu.z * rmu.z);

        rmu.x /= rm;
        rmu.y /= rm;
        rmu.z /= rm;

        rcu.x = coords[wi].x - ctx.topo.solvent_center.x;
        rcu.y = coords[wi].y - ctx.topo.solvent_center.y;
        rcu.z = coords[wi].z - ctx.topo.solvent_center.z;
        rc = sqrt(rcu.x * rcu.x + rcu.y * rcu.y + rcu.z * rcu.z);
        rcu.x /= rc;
        rcu.y /= rc;
        rcu.z /= rc;

        cos_th = rmu.x * rcu.x + rmu.y * rcu.y + rmu.z * rcu.z;
        if (cos_th > 1) {
            cos_th = 1;
        }
        if (cos_th < -1) {
            cos_th = -1;
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
            // ss << wi << ' ' << iis + 1 << '\n';
        }
    }
    // debug(ss.str());
    const real_t PI = 4.0 * atanf(1.0);
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
            ctx.tdum[jmin] = 99999;
        }
    }

    if (iteration != 0 && iteration % itdis_update == 0) {
        for (int is = 0; is < ctx.n_shells; is++) {
            printf("SHELL %d\n", is);
            wshells[is].avtheta /= (float)itdis_update;
            wshells[is].avn_inshell /= (float)itdis_update;
            wshells[is].theta_corr =
                wshells[is].theta_corr + wshells[is].avtheta - acos(wshells[is].cstb);
            printf("average theta = %f, average in shell = %f, theta_corr = %f\n",
                   wshells[is].avtheta * 180 / M_PI, wshells[is].avn_inshell,
                   wshells[is].theta_corr * 180 / M_PI);
            wshells[is].avtheta = 0;
            wshells[is].avn_inshell = 0;
        }
    }

    for (int is = 0; is < ctx.n_shells; is++) {
        if (wshells[is].n_inshell == 0) {
            continue;
        }

        avtdum = 0;
        for (int il = 0; il < wshells[is].n_inshell; il++) {
            ii = ctx.nsort[il][is];
            arg = 1.0f + ((1.0f - 2.0f * (float)(il + 1)) / (float)wshells[is].n_inshell);
            ctx.theta0[il] = acos(arg);
            ctx.theta0[il] = ctx.theta0[il] - 3 * sin(ctx.theta0[il]) * wshells[is].cstb / 2;
            if (ctx.theta0[il] < 0) {
                ctx.theta0[il] = 0;
            }
            if (ctx.theta0[il] > PI) {
                ctx.theta0[il] = PI;
            }



            avtdum += ctx.theta[ii];
            ener = .5 * ctx.md.polarisation_force *
                   pow(ctx.theta[ii] - ctx.theta0[il] + wshells[is].theta_corr, 2);
            ctx.E_restraint.Upolx += ener;

            dv = ctx.md.polarisation_force * (ctx.theta[ii] - ctx.theta0[il] + wshells[is].theta_corr);

            wi = ctx.n_atoms_solute + 3 * ii;

            rmu.x = coords[wi + 1].x + coords[wi + 2].x - 2 * coords[wi].x;
            rmu.y = coords[wi + 1].y + coords[wi + 2].y - 2 * coords[wi].y;
            rmu.z = coords[wi + 1].z + coords[wi + 2].z - 2 * coords[wi].z;

            rm = sqrt(rmu.x * rmu.x + rmu.y * rmu.y + rmu.z * rmu.z);

            rmu.x /= rm;
            rmu.y /= rm;
            rmu.z /= rm;

            rcu.x = coords[wi].x - ctx.topo.solvent_center.x;
            rcu.y = coords[wi].y - ctx.topo.solvent_center.y;
            rcu.z = coords[wi].z - ctx.topo.solvent_center.z;
            rc = sqrt(rcu.x * rcu.x + rcu.y * rcu.y + rcu.z * rcu.z);
            rcu.x /= rc;
            rcu.y /= rc;
            rcu.z /= rc;

            cos_th = rmu.x * rcu.x + rmu.y * rcu.y + rmu.z * rcu.z;
            if (cos_th > 1) {
                cos_th = 1;
            }
            if (cos_th < -1) {
                cos_th = -1;
            }
            f0 = sin(acos(cos_th));
            if (fabs(f0) < (float)k_singular_sin_epsilon) {
                f0 = (float)k_singular_sin_epsilon;
            }
            f0 = -1.0f / f0;
            f0 *= dv;
            ss << ii << ' ' << cos_th << '\n';

            f1O.x = -2 * (rcu.x - rmu.x * cos_th) / rm;
            f1O.y = -2 * (rcu.y - rmu.y * cos_th) / rm;
            f1O.z = -2 * (rcu.z - rmu.z * cos_th) / rm;
            f1H1.x = (rcu.x - rmu.x * cos_th) / rm;
            f1H1.y = (rcu.y - rmu.y * cos_th) / rm;
            f1H1.z = (rcu.z - rmu.z * cos_th) / rm;
            f1H2.x = (rcu.x - rmu.x * cos_th) / rm;
            f1H2.y = (rcu.y - rmu.y * cos_th) / rm;
            f1H2.z = (rcu.z - rmu.z * cos_th) / rm;

            f2.x = (rmu.x - rcu.x * cos_th) / rc;
            f2.y = (rmu.y - rcu.y * cos_th) / rc;
            f2.z = (rmu.z - rcu.z * cos_th) / rc;

            dvelocities[wi].x += f0 * (f1O.x + f2.x);
            dvelocities[wi].y += f0 * (f1O.y + f2.y);
            dvelocities[wi].z += f0 * (f1O.z + f2.z);
            dvelocities[wi + 1].x += f0 * f1H1.x;
            dvelocities[wi + 1].y += f0 * f1H1.y;
            dvelocities[wi + 1].z += f0 * f1H1.z;
            dvelocities[wi + 2].x += f0 * f1H2.x;
            dvelocities[wi + 2].y += f0 * f1H2.y;
            dvelocities[wi + 2].z += f0 * f1H2.z;
        }

        wshells[is].avtheta += avtdum / (float)wshells[is].n_inshell;
        wshells[is].avn_inshell += wshells[is].n_inshell;
    }
    // debug(ss.str());
}
