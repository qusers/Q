#include "system.h"
#include "restraints.h"
#include "utils.h"
#include "math.h"

#include <stdio.h>

void calc_radix_w_forces() {
    double b, db, ener, dv, fexp;
    coord_t dr;
    double shift;

    if (md.radial_force != 0) {
        shift = sqrt(Boltz * Tfree / md.radial_force);
    }
    else {
        shift = 0;
    }

    // Calculate erst and dv. Note all atoms except oxygens are skipped
    for (int i = n_atoms_solute; i < n_atoms; i += 3) {
        dr.x = coords[i].x - topo.solvent_center.x;
        dr.y = coords[i].y - topo.solvent_center.y;
        dr.z = coords[i].z - topo.solvent_center.z;
        b = sqrt(pow(dr.x, 2) + pow(dr.y, 2) + pow(dr.z, 2));
        db = b - (topo.solvent_radius - shift);

        if (db > 0) {
            ener = 0.5 * md.radial_force * pow(db, 2) - Dwmz;
            dv = md.radial_force * db / b;
        }
        else {
            if (b > 0.0) {
                fexp = exp (awmz * db);
                ener = Dwmz * (pow(fexp, 2) - 2 * fexp);
                dv = -2 * Dwmz * awmz * (fexp - pow(fexp, 2)) / b;
            }
            else {
                dv = 0;
                ener = 0;
            }
        }

        // Update energy and forces
        E_restraint.Uradx += ener;
        dvelocities[i].x += dv * dr.x;
        dvelocities[i].y += dv * dr.y;
        dvelocities[i].z += dv * dr.z;
    }

}

void calc_polx_w_forces(int iteration) {
    int wi, imin, jw, ii, iis, jmin;
    double tmin;
    coord_t rmu, rcu, f1O, f1H1, f1H2, f2;
    double rm, rc;
    double cos_th;
    double avtdum, arg, f0, dv;
    double ener;

    for (int is = 0; is < n_shells; is++) {
        wshells[is].n_inshell = 0;
        if (iteration == 0) {
            wshells[is].theta_corr = 0;
        }
    }

    for (int i = 0; i < n_waters; i++) {
        theta[i] = 0;
        theta0[i] = 0;

        wi = n_atoms_solute + 3 * i;

        rmu.x = coords[wi+1].x + coords[wi+2].x - 2 * coords[wi].x;
        rmu.y = coords[wi+1].y + coords[wi+2].y - 2 * coords[wi].y;
        rmu.z = coords[wi+1].z + coords[wi+2].z - 2 * coords[wi].z;

        rm = sqrt(pow(rmu.x, 2) + pow(rmu.y, 2) + pow(rmu.z, 2));

        rmu.x /= rm;
        rmu.y /= rm;
        rmu.z /= rm;

        rcu.x = coords[wi].x - topo.solvent_center.x;
        rcu.y = coords[wi].y - topo.solvent_center.y;
        rcu.z = coords[wi].z - topo.solvent_center.z;
        rc = sqrt(pow(rcu.x, 2) + pow(rcu.y, 2) + pow(rcu.z, 2));
        rcu.x /= rc;
        rcu.y /= rc;
        rcu.z /= rc;

        cos_th = rmu.x * rcu.x + rmu.y*rcu.y + rmu.z*rcu.z;
        if (cos_th > 1) cos_th = 1;
        if (cos_th < -1) cos_th = -1;
        theta[i] = acos(cos_th);
        tdum[i] = theta[i];

        // For waters outside inner shell, locate shell they're in
        if (rc > wshells[n_shells-1].router - wshells[n_shells-1].dr) {
            for (iis = n_shells-1; iis > 0; iis--) {
                if (rc <= wshells[iis].router) break;
            }

            wshells[iis].n_inshell += 1;
            list_sh[wshells[iis].n_inshell-1][iis] = i;
        }
    }

    // Sort the waters according to theta
    for (int is = 0; is < n_shells; is++) {
        imin = 0;
        for (int il = 0; il < wshells[is].n_inshell; il++) {
            tmin = 2 * M_PI;
            for (int jl = 0; jl < wshells[is].n_inshell; jl++) {
                jw = list_sh[jl][is];
                if (tdum[jw] < tmin) {
                    jmin = jw;
                    tmin = theta[jw];
                }
            }
            nsort[imin][is] = jmin;
            imin++;
            tdum[jmin] = 99999;
        }
    }

    // Update theta_corr, averages
    if (iteration != 0 && iteration % itdis_update == 0) {
        for (int is = 0; is < n_shells; is++) {
            printf("SHELL %d\n", is);
            wshells[is].avtheta /= (double) itdis_update;
            wshells[is].avn_inshell /= (double) itdis_update;
            wshells[is].theta_corr = wshells[is].theta_corr + wshells[is].avtheta - acos(wshells[is].cstb);
            printf("average theta = %f, average in shell = %f, theta_corr = %f\n",
                wshells[is].avtheta * 180 / M_PI, wshells[is].avn_inshell, wshells[is].theta_corr * 180 / M_PI);
            wshells[is].avtheta = 0;
            wshells[is].avn_inshell = 0;
        }
    }

    // Calculate energy and force
    for (int is = 0; is < n_shells; is++) {
        if (wshells[is].n_inshell == 0) {
            continue; // Skip empty shell
        }

        avtdum = 0;
        for (int il = 0; il < wshells[is].n_inshell; il++) {
            ii = nsort[il][is];
            arg = 1 + ((1 - 2 * (double) (il+1)) / (double) wshells[is].n_inshell);
            theta0[il] = acos(arg);
            theta0[il] = theta0[il] - 3 * sin(theta0[il]) * wshells[is].cstb / 2;
            if (theta0[il] < 0) theta0[il] = 0;
            if (theta0[il] > M_PI) theta0[il] = M_PI;

            avtdum += theta[ii];
            ener = .5 * md.polarisation_force * pow(theta[ii] - theta0[il] + wshells[is].theta_corr, 2);
            E_restraint.Upolx += ener;

            dv = md.polarisation_force * (theta[ii] - theta0[il] + wshells[is].theta_corr);

            wi = n_atoms_solute + 3 * ii;

            rmu.x = coords[wi+1].x + coords[wi+2].x - 2 * coords[wi].x;
            rmu.y = coords[wi+1].y + coords[wi+2].y - 2 * coords[wi].y;
            rmu.z = coords[wi+1].z + coords[wi+2].z - 2 * coords[wi].z;
    
            rm = sqrt(pow(rmu.x, 2) + pow(rmu.y, 2) + pow(rmu.z, 2));
    
            rmu.x /= rm;
            rmu.y /= rm;
            rmu.z /= rm;
    
            rcu.x = coords[wi].x - topo.solvent_center.x;
            rcu.y = coords[wi].y - topo.solvent_center.y;
            rcu.z = coords[wi].z - topo.solvent_center.z;
            rc = sqrt(pow(rcu.x, 2) + pow(rcu.y, 2) + pow(rcu.z, 2));
            rcu.x /= rc;
            rcu.y /= rc;
            rcu.z /= rc;
    
            cos_th = rmu.x * rcu.x + rmu.y*rcu.y + rmu.z*rcu.z;
            if (cos_th > 1) cos_th = 1;
            if (cos_th < -1) cos_th = -1;
            f0 = sin(acos(cos_th));
            if (abs(f0) < 1.0E-12) f0 = 1.0E-12;
            f0 = -1.0 / f0;
            f0 *= dv;

            f1O.x = -2 * (rcu.x - rmu.x * cos_th) / rm;
            f1O.y = -2 * (rcu.y - rmu.y * cos_th) / rm;
            f1O.z = -2 * (rcu.z - rmu.z * cos_th) / rm;
            f1H1.x =     (rcu.x - rmu.x * cos_th) / rm;
            f1H1.y =     (rcu.y - rmu.y * cos_th) / rm;
            f1H1.z =     (rcu.z - rmu.z * cos_th) / rm;
            f1H2.x =     (rcu.x - rmu.x * cos_th) / rm;
            f1H2.y =     (rcu.y - rmu.y * cos_th) / rm;
            f1H2.z =     (rcu.z - rmu.z * cos_th) / rm;

            f2.x = ( rmu.x - rcu.x * cos_th) / rc;
            f2.y = ( rmu.y - rcu.y * cos_th) / rc;
            f2.z = ( rmu.z - rcu.z * cos_th) / rc;

            dvelocities[wi].x   += f0 * (f1O.x + f2.x);
            dvelocities[wi].y   += f0 * (f1O.y + f2.y);
            dvelocities[wi].z   += f0 * (f1O.z + f2.z);
            dvelocities[wi+1].x += f0 * (f1H1.x);
            dvelocities[wi+1].y += f0 * (f1H1.y);
            dvelocities[wi+1].z += f0 * (f1H1.z);
            dvelocities[wi+2].x += f0 * (f1H2.x);
            dvelocities[wi+2].y += f0 * (f1H2.y);
            dvelocities[wi+2].z += f0 * (f1H2.z);
        }

        wshells[is].avtheta     += avtdum / (double) wshells[is].n_inshell;
        wshells[is].avn_inshell += wshells[is].n_inshell;
    }
}

void calc_pshell_forces() {
    coord_t dr;
    double k, r2, ener;

    for (int i = 0; i < n_atoms_solute; i++) {
        if (shell[i] || excluded[i]) {
            // printf("i = %d excluded = %s shell = %s\n", i, excluded[i] ? "True" : "False", shell[i] ? "True" : "False");
            if (excluded[i]) {
                k = k_fix;
            }
            else {
                k = k_pshell;
            }
            dr.x = coords[i].x - coords_top[i].x;
            dr.y = coords[i].y - coords_top[i].y;
            dr.z = coords[i].z - coords_top[i].z;
            r2 = pow(dr.x, 2) + pow(dr.y, 2) + pow(dr.z, 2);
            ener = 0.5 * k * r2;
            // printf("dr = %f %f %f\n", dr.x, dr.y, dr.z);

            if (excluded[i]) E_restraint.Ufix += ener;
            if (shell[i]) E_restraint.Ushell += ener;

            dvelocities[i].x += k * dr.x;
            dvelocities[i].y += k * dr.y;
            dvelocities[i].z += k * dr.z;
        }
    }
}

// sequence restraints (independent of Q-state)
void calc_restrseq_forces() {
    double k, mass, totmass;
    coord_t dr;
    double r2, ener;

    for (int s = 0; s < n_restrseqs; s++) {
        k = restrseqs[s].k;

        dr.x = 0;
        dr.y = 0;
        dr.z = 0;
        int n_ctr = 0;
        totmass = 0;

        // Geometric center
        if (restrseqs[s].to_center == 1) {
            for (int i = restrseqs[s].ai-1; i < restrseqs[s].aj-1; i++) {
                if (heavy[i] || restrseqs[s].ih) {
                    n_ctr++;
                    dr.x += (coords[i].x - coords_top[i].x);
                    dr.y += (coords[i].y - coords_top[i].y);
                    dr.z += (coords[i].z - coords_top[i].z);
                }
            }

            if (n_ctr > 0) {
                dr.x /= n_ctr;
                dr.y /= n_ctr;
                dr.z /= n_ctr;
                r2 = pow(dr.x, 2) + pow(dr.y, 2) + pow(dr.z, 2);
                ener = .5 * k * r2;
                E_restraint.Upres += ener;

                for (int i = restrseqs[s].ai-1; i < restrseqs[s].aj-1; i++) {
                    if (heavy[i] || restrseqs[s].ih) {
                        mass = catypes[atypes[i].code - 1].m;
                        dvelocities[i].x += (k * dr.x * mass / 12.010);
                        dvelocities[i].y += (k * dr.y * mass / 12.010);
                        dvelocities[i].z += (k * dr.z * mass / 12.010);
                    }
                }
            }
        }

        // Mass center
        else if (restrseqs[s].to_center == 2) {
            for (int i = restrseqs[s].ai-1; i < restrseqs[s].aj-1; i++) {
                if (heavy[i] || restrseqs[i].ih) {
                    mass = catypes[atypes[i].code-1].m;
                    totmass += mass;
                    dr.x += (coords[i].x - coords_top[i].x) * mass;
                    dr.y += (coords[i].y - coords_top[i].y) * mass;
                    dr.z += (coords[i].z - coords_top[i].z) * mass;
                }
            }

            if (totmass > 0) {
                dr.x /= totmass;
                dr.y /= totmass;
                dr.z /= totmass;
                r2 = pow(dr.x, 2) + pow(dr.y, 2) + pow(dr.z, 2);
                ener = .5 * k * r2;
                E_restraint.Upres += ener;

                for (int i = restrseqs[s].ai-1; i < restrseqs[s].aj-1; i++) {
                    if (heavy[i] || restrseqs[s].ih) {
                        mass = catypes[atypes[i].code - 1].m;
                        dvelocities[i].x += k * dr.x;
                        dvelocities[i].y += k * dr.y;
                        dvelocities[i].z += k * dr.z;
                    }
                }
            }
        }

        // Restrain to topology coordinate
        else {
            for (int i = restrseqs[s].ai-1; i < restrseqs[s].aj-1; i++) {
                if (heavy[i] || restrseqs[s].ih) {
                    dr.x = coords[i].x - coords_top[i].x;
                    dr.y = coords[i].y - coords_top[i].y;
                    dr.z = coords[i].z - coords_top[i].z;

                    r2 = pow(dr.x, 2) + pow(dr.y, 2) + pow(dr.z, 2);
                    ener = .5 * k * r2;
                    E_restraint.Upres += ener;

                    dvelocities[i].x += k * dr.x;
                    dvelocities[i].y += k * dr.y;
                    dvelocities[i].z += k * dr.z;
                }
            }
        }
    }
}

// extra positional restraints (Q-state dependent)
void calc_restrpos_forces() {
    int state, i;
    coord_t dr;
    double lambda, ener, x2, y2, z2;

    for (int ir = 0; ir < n_restrspos; ir++) {
        state = restrspos[ir].ipsi-1;
        i = restrspos[ir].a-1;

        dr.x = coords[i].x - restrspos[ir].x.x;
        dr.y = coords[i].y - restrspos[ir].x.y;
        dr.z = coords[i].z - restrspos[ir].x.z;

        if (restrspos[ir].ipsi != 0) {
            lambda = lambdas[state];
        }
        else {
            lambda = 1;
        }

        x2 = pow(dr.x, 2);
        y2 = pow(dr.y, 2);
        z2 = pow(dr.z, 2);

        ener = .5 * restrspos[ir].k.x * x2
            + .5 * restrspos[ir].k.y * y2
            + .5 * restrspos[ir].k.z * z2;

        dvelocities[i].x += restrspos[ir].k.x * dr.x * lambda;
        dvelocities[i].y += restrspos[ir].k.y * dr.y * lambda;
        dvelocities[i].z += restrspos[ir].k.z * dr.z * lambda;

        if (restrspos[ir].ipsi == 0) {
            for (int k = 0; k < n_lambdas; k++) {
                EQ_restraint[k].Urestr += ener;
            }
            if (n_lambdas == 0) {
                E_restraint.Upres += ener;
            }
        }
        else {
            EQ_restraint[state].Urestr += ener;
        }
    }
}

// atom-atom distance restraints (Q-state dependent)
void calc_restrdis_forces() {
    int state, i, j;
    coord_t dr;
    double lambda, b, db, dv, ener;

    for (int ir = 0; ir < n_restrdists; ir++) {
        state =  restrdists[ir].ipsi-1;
        i = restrdists[ir].ai - 1;
        j = restrdists[ir].aj - 1;

        dr.x = coords[j].x - coords[i].x;
        dr.y = coords[j].y - coords[i].y;
        dr.z = coords[j].z - coords[i].z;

        if (restrdists[ir].ipsi != 0) {
            lambda = lambdas[state];
        }
        else {
            lambda = 1;
        }

        b = sqrt(pow(dr.x, 2) + pow(dr.y, 2) + pow(dr.z, 2));
        if (b < restrdists[ir].d1) {
            db = b - restrdists[ir].d1;
        }
        else if (b > restrdists[ir].d2) {
            db = b - restrdists[ir].d2;
        }
        else {
            db = 0;
            continue;
        }

        ener = .5 * restrdists[ir].k * pow(db, 2);
        dv = lambda * restrdists[ir].k * db / b;

        dvelocities[j].x += dr.x * dv;
        dvelocities[j].y += dr.y * dv;
        dvelocities[j].z += dr.z * dv;
        dvelocities[i].x -= dr.x * dv;
        dvelocities[i].y -= dr.y * dv;
        dvelocities[i].z -= dr.z * dv;

        if (restrdists[ir].ipsi == 0) {
            for (int k = 0; k < n_lambdas; k++) {
                EQ_restraint[k].Urestr += ener;
            }
            if (n_lambdas == 0) {
                E_restraint.Upres += ener;
            }
        }
        else {
            EQ_restraint[state].Urestr += ener;
        }
    }
}

// atom-atom-atom angle restraints (Q-state dependent)
void calc_restrang_forces() {
    int state, i, j, k;
    coord_t dr, dr2, di, dk;
    double lambda, r2ij, r2jk, rij, rjk, cos_th, th;
    double dth, dv, ener, f1;

    for (int ir = 0; ir < n_restrangs; ir++) {
        state = restrangs[ir].ipsi-1;
        i = restrangs[ir].ai-1;
        j = restrangs[ir].aj-1;
        k = restrangs[ir].ak-1;

        // distance from atom i to atom j
        dr.x = coords[i].x - coords[j].x;
        dr.y = coords[i].y - coords[j].y;
        dr.z = coords[i].z - coords[j].z;

        // distance from atom k to atom j
        dr.x = coords[k].x - coords[j].x;
        dr.y = coords[k].y - coords[j].y;
        dr.z = coords[k].z - coords[j].z;

        if (restrangs[ir].ipsi != 0) {
            lambda = lambdas[state];
        }
        else {
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
        dth = th - to_radians(restrangs[ir].ang);

        ener = .5 * restrangs[ir].k * pow(dth, 2);
        dv = lambda * restrangs[ir].k * dth;

        f1 = sin(th);
        if (abs(f1) < 1E-12) {
            f1 = -1E-12;
        }
        else {
            f1 = -1 / f1;
        }

        di.x = f1 * (dr2.x / (rij * rjk) - cos_th * dr.x / r2ij );
        di.y = f1 * (dr2.y / (rij * rjk) - cos_th * dr.y / r2ij );
        di.z = f1 * (dr2.z / (rij * rjk) - cos_th * dr.z / r2ij );
        dk.x = f1 * (dr.x / (rij * rjk) - cos_th * dr2.x / r2jk );
        dk.y = f1 * (dr.y / (rij * rjk) - cos_th * dr2.y / r2jk );
        dk.z = f1 * (dr.z / (rij * rjk) - cos_th * dr2.z / r2jk );

        dvelocities[i].x += dv * di.x;
        dvelocities[i].y += dv * di.y;
        dvelocities[i].z += dv * di.z;
        dvelocities[k].x += dv * dk.x;
        dvelocities[k].y += dv * dk.y;
        dvelocities[k].z += dv * dk.z;
        dvelocities[j].x -= dv * (di.x + dk.x);
        dvelocities[j].y -= dv * (di.y + dk.y);
        dvelocities[j].z -= dv * (di.z + dk.z);

        if (restrdists[ir].ipsi == 0) {
            for (int k = 0; k < n_lambdas; k++) {
                EQ_restraint[k].Urestr += ener;
            }
            if (n_lambdas == 0) {
                E_restraint.Upres += ener;
            }
        }
        else {
            EQ_restraint[state].Urestr += ener;
        }
    }
}

// extra half-harmonic wall restraints
void calc_restrwall_forces() {
    double k, b, db, ener, dv, fexp;
    coord_t dr;

    for (int ir = 0; ir < n_restrwalls; ir++) {
        for (int i = restrwalls[ir].ai-1; i < restrwalls[ir].aj-1; i++) {
            if (heavy[i] || restrwalls[ir].ih) {
                dr.x = coords[i].x - topo.solvent_center.x;
                dr.y = coords[i].y - topo.solvent_center.y;
                dr.x = coords[i].x - topo.solvent_center.x;

                b = sqrt(pow(dr.x, 2) + pow(dr.y, 2) + pow(dr.z, 2));
                db = b - restrwalls[ir].d;

                if (db > 0) {
                    ener = .5 * k * pow(db, 2) - restrwalls[ir].dMorse;
                    dv = k * db / b;
                }
                else {
                    fexp = exp(restrwalls[ir].aMorse * db);
                    ener = restrwalls[ir].dMorse * (fexp * fexp - 2 * fexp);
                    dv = -2 * restrwalls[ir].dMorse * restrwalls[ir].aMorse * (fexp - fexp*fexp) / b;
                }
                E_restraint.Upres += ener;

                dvelocities[i].x += dv * dr.x;
                dvelocities[i].y += dv * dr.y;
                dvelocities[i].z += dv * dr.z;
            }
        }
    }
}