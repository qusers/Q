#include "cpu/include/context.h"
#include "cpu/include/constants.h"

#include "system.h"
#include "restraints.h"
#include "utils.h"

#include <math.h>
#include <stdio.h>

void calc_radix_w_forces() {
    auto& ctx = Context::instance();

    double b, db, ener, dv, fexp;
    coord_t dr;
    double shift;

    if (ctx.md.radial_force != 0) {
        shift = sqrt(Boltz * ctx.Tfree / ctx.md.radial_force);
    } else {
        shift = 0;
    }

    // Calculate erst and dv. Note all atoms except oxygens are skipped
    for (int i = ctx.n_atoms_solute; i < ctx.n_atoms; i += 3) {
        dr.x = ctx.coords[i].x - ctx.topo.solvent_center.x;
        dr.y = ctx.coords[i].y - ctx.topo.solvent_center.y;
        dr.z = ctx.coords[i].z - ctx.topo.solvent_center.z;
        b = sqrt(pow(dr.x, 2) + pow(dr.y, 2) + pow(dr.z, 2));
        db = b - (ctx.topo.solvent_radius - shift);

        if (db > 0) {
            ener = 0.5 * ctx.md.radial_force * pow(db, 2) - ctx.Dwmz;
            dv = ctx.md.radial_force * db / b;
        } else {
            if (b > 0.0) {
                fexp = exp(ctx.awmz * db);
                ener = ctx.Dwmz * (pow(fexp, 2) - 2 * fexp);
                dv = -2 * ctx.Dwmz * ctx.awmz * (fexp - pow(fexp, 2)) / b;
            } else {
                dv = 0;
                ener = 0;
            }
        }

        ctx.E_restraint.Uradx += ener;
        ctx.dvelocities[i].x += dv * dr.x;
        ctx.dvelocities[i].y += dv * dr.y;
        ctx.dvelocities[i].z += dv * dr.z;
    }
}

void calc_polx_w_forces(int iteration) {
    auto& ctx = Context::instance();

    int wi, imin, jw, ii, iis, jmin;
    double tmin;
    coord_t rmu, rcu, f1O, f1H1, f1H2, f2;
    double rm, rc;
    double cos_th;
    double avtdum, arg, f0, dv;
    double ener;

    for (int is = 0; is < ctx.n_shells; is++) {
        ctx.wshells[is].n_inshell = 0;
        if (iteration == 0) {
            ctx.wshells[is].theta_corr = 0;
        }
    }

    for (int i = 0; i < ctx.n_waters; i++) {
        ctx.theta[i] = 0;
        ctx.theta0[i] = 0;

        wi = ctx.n_atoms_solute + 3 * i;

        rmu.x = ctx.coords[wi + 1].x + ctx.coords[wi + 2].x - 2 * ctx.coords[wi].x;
        rmu.y = ctx.coords[wi + 1].y + ctx.coords[wi + 2].y - 2 * ctx.coords[wi].y;
        rmu.z = ctx.coords[wi + 1].z + ctx.coords[wi + 2].z - 2 * ctx.coords[wi].z;

        rm = sqrt(pow(rmu.x, 2) + pow(rmu.y, 2) + pow(rmu.z, 2));

        rmu.x /= rm;
        rmu.y /= rm;
        rmu.z /= rm;

        rcu.x = ctx.coords[wi].x - ctx.topo.solvent_center.x;
        rcu.y = ctx.coords[wi].y - ctx.topo.solvent_center.y;
        rcu.z = ctx.coords[wi].z - ctx.topo.solvent_center.z;
        rc = sqrt(pow(rcu.x, 2) + pow(rcu.y, 2) + pow(rcu.z, 2));
        rcu.x /= rc;
        rcu.y /= rc;
        rcu.z /= rc;

        cos_th = rmu.x * rcu.x + rmu.y * rcu.y + rmu.z * rcu.z;
        if (cos_th > 1) cos_th = 1;
        if (cos_th < -1) cos_th = -1;
        ctx.theta[i] = acos(cos_th);
        ctx.tdum[i] = ctx.theta[i];

        if (rc > ctx.wshells[ctx.n_shells - 1].router - ctx.wshells[ctx.n_shells - 1].dr) {
            for (iis = ctx.n_shells - 1; iis > 0; iis--) {
                if (rc <= ctx.wshells[iis].router) break;
            }

            ctx.wshells[iis].n_inshell += 1;
            ctx.list_sh[ctx.wshells[iis].n_inshell - 1][iis] = i;
        }
    }

    for (int is = 0; is < ctx.n_shells; is++) {
        imin = 0;
        for (int il = 0; il < ctx.wshells[is].n_inshell; il++) {
            tmin = 2 * M_PI;
            for (int jl = 0; jl < ctx.wshells[is].n_inshell; jl++) {
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
            ctx.wshells[is].avtheta /= (double)itdis_update;
            ctx.wshells[is].avn_inshell /= (double)itdis_update;
            ctx.wshells[is].theta_corr =
                ctx.wshells[is].theta_corr + ctx.wshells[is].avtheta - acos(ctx.wshells[is].cstb);
            printf("average theta = %f, average in shell = %f, theta_corr = %f\n",
                   ctx.wshells[is].avtheta * 180 / M_PI,
                   ctx.wshells[is].avn_inshell,
                   ctx.wshells[is].theta_corr * 180 / M_PI);
            ctx.wshells[is].avtheta = 0;
            ctx.wshells[is].avn_inshell = 0;
        }
    }

    for (int is = 0; is < ctx.n_shells; is++) {
        if (ctx.wshells[is].n_inshell == 0) {
            continue;
        }

        avtdum = 0;
        for (int il = 0; il < ctx.wshells[is].n_inshell; il++) {
            ii = ctx.nsort[il][is];
            arg = 1 + ((1 - 2 * (double)(il + 1)) / (double)ctx.wshells[is].n_inshell);
            ctx.theta0[il] = acos(arg);
            ctx.theta0[il] = ctx.theta0[il] - 3 * sin(ctx.theta0[il]) * ctx.wshells[is].cstb / 2;
            if (ctx.theta0[il] < 0) ctx.theta0[il] = 0;
            if (ctx.theta0[il] > M_PI) ctx.theta0[il] = M_PI;

            avtdum += ctx.theta[ii];
            ener = .5 * ctx.md.polarisation_force *
                   pow(ctx.theta[ii] - ctx.theta0[il] + ctx.wshells[is].theta_corr, 2);
            ctx.E_restraint.Upolx += ener;

            dv = ctx.md.polarisation_force * (ctx.theta[ii] - ctx.theta0[il] + ctx.wshells[is].theta_corr);

            wi = ctx.n_atoms_solute + 3 * ii;

            rmu.x = ctx.coords[wi + 1].x + ctx.coords[wi + 2].x - 2 * ctx.coords[wi].x;
            rmu.y = ctx.coords[wi + 1].y + ctx.coords[wi + 2].y - 2 * ctx.coords[wi].y;
            rmu.z = ctx.coords[wi + 1].z + ctx.coords[wi + 2].z - 2 * ctx.coords[wi].z;

            rm = sqrt(pow(rmu.x, 2) + pow(rmu.y, 2) + pow(rmu.z, 2));

            rmu.x /= rm;
            rmu.y /= rm;
            rmu.z /= rm;

            rcu.x = ctx.coords[wi].x - ctx.topo.solvent_center.x;
            rcu.y = ctx.coords[wi].y - ctx.topo.solvent_center.y;
            rcu.z = ctx.coords[wi].z - ctx.topo.solvent_center.z;
            rc = sqrt(pow(rcu.x, 2) + pow(rcu.y, 2) + pow(rcu.z, 2));
            rcu.x /= rc;
            rcu.y /= rc;
            rcu.z /= rc;

            cos_th = rmu.x * rcu.x + rmu.y * rcu.y + rmu.z * rcu.z;
            if (cos_th > 1) cos_th = 1;
            if (cos_th < -1) cos_th = -1;
            f0 = sin(acos(cos_th));
            if (fabs(f0) < 1.0E-12) f0 = 1.0E-12;
            f0 = -1.0 / f0;
            f0 *= dv;

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

            ctx.dvelocities[wi].x += f0 * (f1O.x + f2.x);
            ctx.dvelocities[wi].y += f0 * (f1O.y + f2.y);
            ctx.dvelocities[wi].z += f0 * (f1O.z + f2.z);
            ctx.dvelocities[wi + 1].x += f0 * f1H1.x;
            ctx.dvelocities[wi + 1].y += f0 * f1H1.y;
            ctx.dvelocities[wi + 1].z += f0 * f1H1.z;
            ctx.dvelocities[wi + 2].x += f0 * f1H2.x;
            ctx.dvelocities[wi + 2].y += f0 * f1H2.y;
            ctx.dvelocities[wi + 2].z += f0 * f1H2.z;
        }

        ctx.wshells[is].avtheta += avtdum / (double)ctx.wshells[is].n_inshell;
        ctx.wshells[is].avn_inshell += ctx.wshells[is].n_inshell;
    }
}

void calc_pshell_forces() {
    auto& ctx = Context::instance();

    coord_t dr;
    double k, r2, ener;

    for (int i = 0; i < ctx.n_atoms_solute; i++) {
        if (ctx.shell[i] || ctx.excluded[i]) {
            if (ctx.excluded[i]) {
                k = k_fix;
            } else {
                k = k_pshell;
            }

            dr.x = ctx.coords[i].x - ctx.coords_top[i].x;
            dr.y = ctx.coords[i].y - ctx.coords_top[i].y;
            dr.z = ctx.coords[i].z - ctx.coords_top[i].z;
            r2 = pow(dr.x, 2) + pow(dr.y, 2) + pow(dr.z, 2);
            ener = 0.5 * k * r2;

            if (ctx.excluded[i]) ctx.E_restraint.Ufix += ener;
            if (ctx.shell[i]) ctx.E_restraint.Ushell += ener;

            ctx.dvelocities[i].x += k * dr.x;
            ctx.dvelocities[i].y += k * dr.y;
            ctx.dvelocities[i].z += k * dr.z;
        }
    }
}

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
                    dr.x += (ctx.coords[i].x - ctx.coords_top[i].x);
                    dr.y += (ctx.coords[i].y - ctx.coords_top[i].y);
                    dr.z += (ctx.coords[i].z - ctx.coords_top[i].z);
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
                    dr.x += (ctx.coords[i].x - ctx.coords_top[i].x) * mass;
                    dr.y += (ctx.coords[i].y - ctx.coords_top[i].y) * mass;
                    dr.z += (ctx.coords[i].z - ctx.coords_top[i].z) * mass;
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
                    dr.x = ctx.coords[i].x - ctx.coords_top[i].x;
                    dr.y = ctx.coords[i].y - ctx.coords_top[i].y;
                    dr.z = ctx.coords[i].z - ctx.coords_top[i].z;

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

void calc_restrpos_forces() {
    auto& ctx = Context::instance();

    int state, i;
    coord_t dr;
    double lambda, ener, x2, y2, z2;

    for (int ir = 0; ir < ctx.n_restrspos; ir++) {
        state = ctx.restrspos[ir].ipsi - 1;
        i = ctx.restrspos[ir].a - 1;

        dr.x = ctx.coords[i].x - ctx.restrspos[ir].x.x;
        dr.y = ctx.coords[i].y - ctx.restrspos[ir].x.y;
        dr.z = ctx.coords[i].z - ctx.restrspos[ir].x.z;

        if (ctx.restrspos[ir].ipsi != 0) {
            lambda = ctx.lambdas[state];
        } else {
            lambda = 1;
        }

        x2 = pow(dr.x, 2);
        y2 = pow(dr.y, 2);
        z2 = pow(dr.z, 2);

        ener = .5 * ctx.restrspos[ir].k.x * x2 + .5 * ctx.restrspos[ir].k.y * y2 +
               .5 * ctx.restrspos[ir].k.z * z2;

        ctx.dvelocities[i].x += ctx.restrspos[ir].k.x * dr.x * lambda;
        ctx.dvelocities[i].y += ctx.restrspos[ir].k.y * dr.y * lambda;
        ctx.dvelocities[i].z += ctx.restrspos[ir].k.z * dr.z * lambda;

        if (ctx.restrspos[ir].ipsi == 0) {
            for (int k = 0; k < ctx.n_lambdas; k++) {
                ctx.EQ_restraint[k].Urestr += ener;
            }
            if (ctx.n_lambdas == 0) {
                ctx.E_restraint.Upres += ener;
            }
        } else {
            ctx.EQ_restraint[state].Urestr += ener;
        }
    }
}

void calc_restrdis_forces() {
    auto& ctx = Context::instance();

    int state, i, j;
    coord_t dr;
    double lambda, b, db, dv, ener;

    for (int ir = 0; ir < ctx.n_restrdists; ir++) {
        state = ctx.restrdists[ir].ipsi - 1;
        i = ctx.restrdists[ir].ai - 1;
        j = ctx.restrdists[ir].aj - 1;

        dr.x = ctx.coords[j].x - ctx.coords[i].x;
        dr.y = ctx.coords[j].y - ctx.coords[i].y;
        dr.z = ctx.coords[j].z - ctx.coords[i].z;

        if (ctx.restrdists[ir].ipsi != 0) {
            lambda = ctx.lambdas[state];
        } else {
            lambda = 1;
        }

        b = sqrt(pow(dr.x, 2) + pow(dr.y, 2) + pow(dr.z, 2));
        if (b < ctx.restrdists[ir].d1) {
            db = b - ctx.restrdists[ir].d1;
        } else if (b > ctx.restrdists[ir].d2) {
            db = b - ctx.restrdists[ir].d2;
        } else {
            db = 0;
            continue;
        }

        ener = .5 * ctx.restrdists[ir].k * pow(db, 2);
        dv = lambda * ctx.restrdists[ir].k * db / b;

        ctx.dvelocities[j].x += dr.x * dv;
        ctx.dvelocities[j].y += dr.y * dv;
        ctx.dvelocities[j].z += dr.z * dv;
        ctx.dvelocities[i].x -= dr.x * dv;
        ctx.dvelocities[i].y -= dr.y * dv;
        ctx.dvelocities[i].z -= dr.z * dv;

        if (ctx.restrdists[ir].ipsi == 0) {
            for (int k = 0; k < ctx.n_lambdas; k++) {
                ctx.EQ_restraint[k].Urestr += ener;
            }
            if (ctx.n_lambdas == 0) {
                ctx.E_restraint.Upres += ener;
            }
        } else {
            ctx.EQ_restraint[state].Urestr += ener;
        }
    }
}

void calc_restrang_forces() {
    auto& ctx = Context::instance();

    int state, i, j, k;
    coord_t dr, dr2, di, dk;
    double lambda, r2ij, r2jk, rij, rjk, cos_th, th;
    double dth, dv, ener, f1;

    for (int ir = 0; ir < ctx.n_restrangs; ir++) {
        state = ctx.restrangs[ir].ipsi - 1;
        i = ctx.restrangs[ir].ai - 1;
        j = ctx.restrangs[ir].aj - 1;
        k = ctx.restrangs[ir].ak - 1;

        dr.x = ctx.coords[i].x - ctx.coords[j].x;
        dr.y = ctx.coords[i].y - ctx.coords[j].y;
        dr.z = ctx.coords[i].z - ctx.coords[j].z;

        dr2.x = ctx.coords[k].x - ctx.coords[j].x;
        dr2.y = ctx.coords[k].y - ctx.coords[j].y;
        dr2.z = ctx.coords[k].z - ctx.coords[j].z;

        if (ctx.restrangs[ir].ipsi != 0) {
            lambda = ctx.lambdas[state];
        } else {
            lambda = 1;
        }

        r2ij = pow(dr.x, 2) + pow(dr.y, 2) + pow(dr.z, 2);
        r2jk = pow(dr2.x, 2) + pow(dr2.y, 2) + pow(dr2.z, 2);

        rij = sqrt(r2ij);
        rjk = sqrt(r2jk);

        cos_th = dr.x * dr2.x + dr.y * dr2.y + dr.z * dr2.z;
        cos_th /= rij * rjk;

        if (cos_th > 1) cos_th = 1;
        if (cos_th < -1) cos_th = -1;

        th = acos(cos_th);
        dth = th - to_radians(ctx.restrangs[ir].ang);

        ener = .5 * ctx.restrangs[ir].k * pow(dth, 2);
        dv = lambda * ctx.restrangs[ir].k * dth;

        f1 = sin(th);
        if (fabs(f1) < 1E-12) {
            f1 = -1E-12;
        } else {
            f1 = -1 / f1;
        }

        di.x = f1 * (dr2.x / (rij * rjk) - cos_th * dr.x / r2ij);
        di.y = f1 * (dr2.y / (rij * rjk) - cos_th * dr.y / r2ij);
        di.z = f1 * (dr2.z / (rij * rjk) - cos_th * dr.z / r2ij);
        dk.x = f1 * (dr.x / (rij * rjk) - cos_th * dr2.x / r2jk);
        dk.y = f1 * (dr.y / (rij * rjk) - cos_th * dr2.y / r2jk);
        dk.z = f1 * (dr.z / (rij * rjk) - cos_th * dr2.z / r2jk);

        ctx.dvelocities[i].x += dv * di.x;
        ctx.dvelocities[i].y += dv * di.y;
        ctx.dvelocities[i].z += dv * di.z;
        ctx.dvelocities[k].x += dv * dk.x;
        ctx.dvelocities[k].y += dv * dk.y;
        ctx.dvelocities[k].z += dv * dk.z;
        ctx.dvelocities[j].x -= dv * (di.x + dk.x);
        ctx.dvelocities[j].y -= dv * (di.y + dk.y);
        ctx.dvelocities[j].z -= dv * (di.z + dk.z);

        if (ctx.restrangs[ir].ipsi == 0) {
            for (int lambda_idx = 0; lambda_idx < ctx.n_lambdas; lambda_idx++) {
                ctx.EQ_restraint[lambda_idx].Urestr += ener;
            }
            if (ctx.n_lambdas == 0) {
                ctx.E_restraint.Upres += ener;
            }
        } else {
            ctx.EQ_restraint[state].Urestr += ener;
        }
    }
}

void calc_restrwall_forces() {
    auto& ctx = Context::instance();

    double k, b, db, ener, dv, fexp;
    coord_t dr;

    for (int ir = 0; ir < ctx.n_restrwalls; ir++) {
        k = ctx.restrwalls[ir].k;
        for (int i = ctx.restrwalls[ir].ai - 1; i < ctx.restrwalls[ir].aj - 1; i++) {
            if (ctx.heavy[i] || ctx.restrwalls[ir].ih) {
                dr.x = ctx.coords[i].x - ctx.topo.solvent_center.x;
                dr.y = ctx.coords[i].y - ctx.topo.solvent_center.y;
                dr.z = ctx.coords[i].z - ctx.topo.solvent_center.z;

                b = sqrt(pow(dr.x, 2) + pow(dr.y, 2) + pow(dr.z, 2));
                db = b - ctx.restrwalls[ir].d;

                if (db > 0) {
                    ener = .5 * k * pow(db, 2) - ctx.restrwalls[ir].dMorse;
                    dv = k * db / b;
                } else {
                    fexp = exp(ctx.restrwalls[ir].aMorse * db);
                    ener = ctx.restrwalls[ir].dMorse * (fexp * fexp - 2 * fexp);
                    dv = -2 * ctx.restrwalls[ir].dMorse * ctx.restrwalls[ir].aMorse *
                         (fexp - fexp * fexp) / b;
                }
                ctx.E_restraint.Upres += ener;

                ctx.dvelocities[i].x += dv * dr.x;
                ctx.dvelocities[i].y += dv * dr.y;
                ctx.dvelocities[i].z += dv * dr.z;
            }
        }
    }
}
