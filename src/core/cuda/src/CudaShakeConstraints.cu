#include <iostream>

#include "cuda/include/CudaContext.cuh"
#include "cuda/include/CudaShakeConstraints.cuh"
#include "utils.h"

__global__ void calc_shake_constraints_kernel(
    int n_molecules,
    int* mol_n_shakes,
    shake_bond_t* shake_bonds,
    coord_t* coords,
    coord_t* xcoords,
    double* winv,
    int* total_iterations,
    int* mol_shake_offset) {
    int idx = blockIdx.x;
    if (idx >= n_molecules) return;

    int mol = idx;

    int ai, aj, n_iterations, shake;
    double xij2, diff, corr, scp, xxij2;
    coord_t xij, xxij;

    if (mol_n_shakes[mol] == 0) return;
    shake = mol_shake_offset[mol];
    n_iterations = 0;

    bool converged = false;
    if (threadIdx.x == 0) {
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
            return;
        }

        atomicAdd(total_iterations, n_iterations);
    }
}

void init_shake_constraints_data() {
    int* mol_shake_offset_host = (int*)malloc(sizeof(int) * n_molecules);
    mol_shake_offset_host[0] = 0;
    for (int i = 1; i < n_molecules; i++) {
        mol_shake_offset_host[i] = mol_shake_offset_host[i - 1] + mol_n_shakes[i - 1];
    }

    CudaContext& ctx = CudaContext::instance();
    cudaMemcpy(ctx.d_mol_shake_offset, mol_shake_offset_host, sizeof(int) * n_molecules, cudaMemcpyHostToDevice);
    free(mol_shake_offset_host);
}

int calc_shake_constraints_host() {
    int* d_total_iterations;
    check_cudaMalloc((void**)&d_total_iterations, sizeof(int));
    int total_iterations_host = 0;
    cudaMemcpy(d_total_iterations, &total_iterations_host, sizeof(int), cudaMemcpyHostToDevice);

    int blocks = n_molecules;
    int threads = 32;

    CudaContext& ctx = CudaContext::instance();
    auto d_mol_n_shakes = ctx.d_mol_n_shakes;
    auto d_shake_bonds = ctx.d_shake_bonds;
    auto d_coords = ctx.d_coords;
    auto d_xcoords = ctx.d_xcoords;
    auto d_winv = ctx.d_winv;
    auto d_mol_shake_offset = ctx.d_mol_shake_offset;

    calc_shake_constraints_kernel<<<blocks, threads>>>(
        n_molecules,
        d_mol_n_shakes,
        d_shake_bonds,
        d_coords,
        d_xcoords,
        d_winv,
        d_total_iterations,
        d_mol_shake_offset);
    cudaDeviceSynchronize();
    cudaMemcpy(&total_iterations_host, d_total_iterations, sizeof(int), cudaMemcpyDeviceToHost);
    cudaMemcpy(coords, d_coords, sizeof(coord_t) * n_atoms, cudaMemcpyDeviceToHost);
    return total_iterations_host;
}
