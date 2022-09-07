#include "system.h"
#include "shake.h"

#include <math.h>
#include <stdio.h>

int calc_shake_constraints(coord_t *coords, coord_t *xcoords) {
    int ai, aj, n_iterations, total_iterations, shake;
    double xij2, diff, corr, scp, xxij2;
    coord_t xij, xxij;

    shake = 0;
    for (int mol = 0; mol < n_molecules; mol++) {
        if (mol_n_shakes[mol] == 0) continue;
        n_iterations = 0;

        bool converged = false;
        do {
            for (int i = 0; i < mol_n_shakes[mol]; i++) {
                shake_bonds[shake + i].ready = false;
            }
    
            for (int i = 0; i < mol_n_shakes[mol]; i++) {
                if (!shake_bonds[shake + i].ready) {
                    ai = shake_bonds[shake + i].ai-1;
                    aj = shake_bonds[shake + i].aj-1;
    
                    xij.x = coords[ai].x - coords[aj].x;
                    xij.y = coords[ai].y - coords[aj].y;
                    xij.z = coords[ai].z - coords[aj].z;
                    xij2 = pow(xij.x, 2) + pow(xij.y, 2) + pow(xij.z, 2);
                    diff = shake_bonds[shake + i].dist2 - xij2;
                    if (abs(diff) < shake_tol * shake_bonds[shake + i].dist2) {
                        shake_bonds[shake + i].ready = true;
                    }
                    xxij.x = xcoords[ai].x - xcoords[aj].x;
                    xxij.y = xcoords[ai].y - xcoords[aj].y;
                    xxij.z = xcoords[ai].z - xcoords[aj].z;
                    scp = xij.x * xxij.x + xij.y * xxij.y + xij.z * xxij.z;
                    corr = diff / (2 * scp * (winv[ai] + winv[aj]));
    
                    coords[ai].x += xxij.x * corr * winv[ai];
                    coords[ai].y += xxij.y * corr * winv[ai];
                    coords[ai].z += xxij.z * corr * winv[ai];
                    coords[aj].x -= xxij.x * corr * winv[aj];
                    coords[aj].y -= xxij.y * corr * winv[aj];
                    coords[aj].z -= xxij.z * corr * winv[aj];
                }
            }
    
            n_iterations++;
    
            converged = true;
            for (int i = 0; i < mol_n_shakes[mol]; i++) {
                if (!shake_bonds[shake + i].ready) {
                    converged = false;
                    break;
                }
            }
        } while(n_iterations < shake_max_iter && !converged);

        if (!converged) {
            for (int i = 0; i < mol_n_shakes[mol]; i++) {
                ai = shake_bonds[shake + i].ai-1;
                aj = shake_bonds[shake + i].aj-1;

                xxij.x = xcoords[ai].x - xcoords[aj].x;
                xxij.y = xcoords[ai].y - xcoords[aj].y;
                xxij.z = xcoords[ai].z - xcoords[aj].z;
                xxij2 = pow(xxij.x, 2) + pow(xxij.y, 2) + pow(xxij.z, 2);
                printf(">>> Shake failed, i = %d,j = %d, d = %f, d0 = %f", ai, aj, sqrt(xxij2), shake_bonds[shake + i].dist2);
            }
            exit(EXIT_FAILURE);
        }

        shake += mol_n_shakes[mol];
        total_iterations += n_iterations;
    }

    // Set niter to the average number of iterations per molecule
    return total_iterations / n_molecules;
}

void initial_shaking() {
    for (int i = 0; i < n_atoms; i++) {
        xcoords[i].x = coords[i].x;
        xcoords[i].y = coords[i].y;
        xcoords[i].z = coords[i].z;
    }
    calc_shake_constraints(coords, xcoords);
    for (int i = 0; i < n_atoms; i++) {
        xcoords[i].x = coords[i].x - dt * velocities[i].x;
        xcoords[i].y = coords[i].y - dt * velocities[i].y;
        xcoords[i].z = coords[i].z - dt * velocities[i].z;
    }
    calc_shake_constraints(xcoords, coords);
    for (int i = 0; i < n_atoms; i++) {
        velocities[i].x = (coords[i].x - xcoords[i].x) / dt;
        velocities[i].y = (coords[i].y - xcoords[i].y) / dt;
        velocities[i].z = (coords[i].z - xcoords[i].z) / dt;
    }
}

void stop_cm_translation() {
    double total_mass = 0;
    double rmass = 0;
    coord_t vcm;
    vcm.x = 0;
    vcm.y = 0;
    vcm.z = 0;

    for (int ai = 0; ai < n_atoms; ai++) {
        rmass = catypes[atypes[ai].code-1].m;
        total_mass += rmass;
        vcm.x = vcm.x + velocities[ai].x * rmass;
        vcm.y = vcm.y + velocities[ai].y;
        vcm.z = vcm.z + velocities[ai].z;
    }

    vcm.x = vcm.x / total_mass;
    vcm.y = vcm.y / total_mass;
    vcm.z = vcm.z / total_mass;

    for (int ai = 0; ai < n_atoms; ai++) {
        velocities[ai].x = velocities[ai].x - vcm.x;
        velocities[ai].y = velocities[ai].y - vcm.y;
        velocities[ai].z = velocities[ai].z - vcm.z;
    }
}