#include "cuda/include/CudaShakeConstraints.cuh"

__global__ void calc_shake_constraints_kernel() {
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
                    ai = shake_bonds[shake + i].ai - 1;
                    aj = shake_bonds[shake + i].aj - 1;

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
        } while (n_iterations < shake_max_iter && !converged);

        if (!converged) {
            for (int i = 0; i < mol_n_shakes[mol]; i++) {
                ai = shake_bonds[shake + i].ai - 1;
                aj = shake_bonds[shake + i].aj - 1;

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