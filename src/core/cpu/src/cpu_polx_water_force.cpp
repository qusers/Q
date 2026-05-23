#include "cpu_polx_water_force.h"

#include <math.h>
#include <stdio.h>
#include <stdlib.h>

#include "constants.h"
#include "context.h"

void calc_polx_w_forces(int iteration) {
    auto& ctx = Context::instance();
    auto &coords = ctx.coords->cpu_data_p;
    auto &dvelocities = ctx.dvelocities->cpu_data_p;
    auto *excluded = ctx.excluded->cpu_data_p;
    auto *wshells = ctx.wshells->cpu_data_p;

    int wi, imin, jw, ii, iis, jmin;
    real_t tmin;
    coord_t rmu, rcu, f1O, f1H1, f1H2, f2;
    real_t rm, rc;
    real_t cos_th;
    real_t avtdum, arg, f0, dv;
    real_t ener;

    for (int is = 0; is < ctx.n_shells; is++) {
        wshells[is].n_inshell = 0;
        if (iteration == 0 && !ctx.has_restart_wshell_theta_corr) {
            wshells[is].theta_corr = 0;
        }
    }

    for (int i = 0; i < ctx.n_waters; i++) {
        ctx.theta[i] = 0;
        ctx.theta0[i] = 0;

        wi = ctx.n_atoms_solute + 3 * i;
        if (excluded[wi]) continue;

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
        }
    }

    for (int is = 0; is < ctx.n_shells; is++) {
        imin = 0;
        for (int il = 0; il < wshells[is].n_inshell; il++) {
            tmin = 2 * static_cast<real_t>(q_fortran_pi);
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
            const float avtheta_f = static_cast<float>(wshells[is].avtheta) / static_cast<float>(itdis_update);
            wshells[is].avtheta = static_cast<real_t>(avtheta_f);
            wshells[is].avn_inshell /= static_cast<real_t>(itdis_update);
            wshells[is].theta_corr =
                static_cast<real_t>(static_cast<float>(
                    static_cast<float>(wshells[is].theta_corr) + avtheta_f -
                    acosf(static_cast<float>(wshells[is].cstb))));
            if (const char* debug = getenv("QGPU_WPOL_DEBUG_UPDATE")) {
                if (debug[0] != '\0') {
                    printf("QGPU_WPOL_UPDATE iter=%d shell=%d avtheta=%.17g avn_inshell=%.17g cstb=%.17g acos_cstb=%.17g theta_corr=%.17g\n",
                           iteration,
                           is + 1,
                           static_cast<double>(wshells[is].avtheta),
                           static_cast<double>(wshells[is].avn_inshell),
                           static_cast<double>(wshells[is].cstb),
                           static_cast<double>(static_cast<real_t>(acosf(static_cast<float>(wshells[is].cstb)))),
                           static_cast<double>(wshells[is].theta_corr));
                }
            }
            printf("average theta = %f, average in shell = %f, theta_corr = %f\n",
                   wshells[is].avtheta * 180 / static_cast<real_t>(q_fortran_pi), wshells[is].avn_inshell,
                   wshells[is].theta_corr * 180 / static_cast<real_t>(q_fortran_pi));
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
            const float arg_f = 1.0f + ((1.0f - 2.0f * static_cast<float>(il + 1)) /
                                        static_cast<float>(wshells[is].n_inshell));
            arg = static_cast<real_t>(arg_f);
            ctx.theta0[il] = acos(arg);
            ctx.theta0[il] = ctx.theta0[il] - 3 * sin(ctx.theta0[il]) * wshells[is].cstb / 2;
            if (ctx.theta0[il] < 0) {
                ctx.theta0[il] = 0;
            }
            if (ctx.theta0[il] > static_cast<real_t>(q_fortran_pi)) {
                ctx.theta0[il] = static_cast<real_t>(q_fortran_pi);
            }

            avtdum += ctx.theta[ii];
            const real_t dtheta = ctx.theta[ii] - ctx.theta0[il] + wshells[is].theta_corr;
            ener = static_cast<real_t>(0.5) * ctx.md.polarisation_force * dtheta * dtheta;
            ctx.E_restraint.Upolx += ener;

            dv = ctx.md.polarisation_force * dtheta;

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
            if (fabs(f0) < real_t(1.0e-12)) {
                f0 = real_t(1.0e-12);
            }
            f0 = static_cast<real_t>(-1.0) / f0;
            f0 = dv * f0;

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

        const real_t shell_avg_theta = avtdum / static_cast<real_t>(wshells[is].n_inshell);
        if (const char* debug = getenv("QGPU_WPOL_DEBUG_AVG")) {
            if (debug[0] != '\0') {
                printf("QGPU_WPOL_AVG iter=%d shell=%d n=%d avtdum=%.17g avg=%.17g avtheta_before=%.17g\n",
                       iteration,
                       is + 1,
                       wshells[is].n_inshell,
                       static_cast<double>(avtdum),
                       static_cast<double>(shell_avg_theta),
                       static_cast<double>(wshells[is].avtheta));
            }
        }
        wshells[is].avtheta = static_cast<real_t>(static_cast<float>(
            static_cast<float>(wshells[is].avtheta) + shell_avg_theta));
        wshells[is].avn_inshell += wshells[is].n_inshell;
    }
}
