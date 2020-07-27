#include "system.h"
#include "restraints.h"

#include <stdio.h>

void calc_radix_w_forces() {
    double b, db, ener, dv, fexp;
    coord_t dr;
    double shift;

    if (k_wsphere != 0) {
        shift = sqrt(Boltz * Temp / k_wsphere);
    }
    else {
        shift = 0;
    }

    // Initialize radial restriction energy
    energies.Uradx = 0;
    
    // Calculate erst and dv. Note all atoms except oxygens are skipped
    for (int i = n_atoms_solute; i < n_atoms; i += 3) {
        dr.x = coords[i].x - centerX;
        dr.y = coords[i].y - centerY;
        dr.z = coords[i].z - centerZ;
        b = sqrt(pow(dr.x, 2) + pow(dr.y, 2) + pow(dr.z, 2));
        db = b - (rwater - shift);

        if (db > 0) {
            ener = 0.5 * k_wsphere * pow(db, 2) - Dwmz;
            dv = k_wsphere * db / b;
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
        energies.Uradx += ener;
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

    // Initialize polar restriction energy
    energies.Upolx = 0;
    
    for (int is = 0; is < n_shells; is++) {
        wshells[is].n_inshell = 0;
        wshells[is].theta_corr = 0;
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

        rcu.x = coords[wi].x - centerX;
        rcu.y = coords[wi].y - centerY;
        rcu.z = coords[wi].z - centerZ;
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
            ener = .5 * k_wpol * pow(theta[ii] - theta0[il] + wshells[is].theta_corr, 2);
            energies.Upolx += ener;

            dv = k_wpol * (theta[ii] - theta0[il] + wshells[is].theta_corr);

            wi = n_atoms_solute + 3 * ii;

            rmu.x = coords[wi+1].x + coords[wi+2].x - 2 * coords[wi].x;
            rmu.y = coords[wi+1].y + coords[wi+2].y - 2 * coords[wi].y;
            rmu.z = coords[wi+1].z + coords[wi+2].z - 2 * coords[wi].z;
    
            rm = sqrt(pow(rmu.x, 2) + pow(rmu.y, 2) + pow(rmu.z, 2));
    
            rmu.x /= rm;
            rmu.y /= rm;
            rmu.z /= rm;
    
            rcu.x = coords[wi].x - centerX;
            rcu.y = coords[wi].y - centerY;
            rcu.z = coords[wi].z - centerZ;
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

    energies.Ushell = 0;
    energies.Ufix = 0;
    for (int i = 0; i < n_atoms_solute; i++) {
        if (excluded[i] || shell[i]) {
            if (excluded[i]) {
                k = k_fix;
            }
            else {
                printf("i = %d\n", i);
                k = k_pshell;
            }
            dr.x = coords[i].x - coords_top[i].x;
            dr.y = coords[i].y - coords_top[i].y;
            dr.z = coords[i].z - coords_top[i].z;
            r2 = pow(dr.x, 2) + pow(dr.y, 2) + pow(dr.z, 2);
            ener = 0.5 * k * r2;

            if (excluded[i]) energies.Ufix += ener;
            if (shell[i]) energies.Ushell += ener;

            dvelocities[i].x += k * dr.x;
            dvelocities[i].y += k * dr.y;
            dvelocities[i].z += k * dr.z;
        }
    }
}

void calc_restrseq_forces() {
    double k, mass, totmass;
    coord_t dr;
    double r2, ener;

    energies.Upres = 0;
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
                energies.Upres += ener;

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
                energies.Upres += ener;

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
                    energies.Upres += ener;

                    dvelocities[i].x += k * dr.x;
                    dvelocities[i].y += k * dr.y;
                    dvelocities[i].z += k * dr.z;
                }
            }
        }
    }
}