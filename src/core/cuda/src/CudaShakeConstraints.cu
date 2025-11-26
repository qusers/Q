#include <iostream>

#include "cuda/include/CudaShakeConstraints.cuh"
#include "utils.h"

namespace CudaShakeConstraints {
bool is_initialized = false;
int* d_mol_n_shakes = nullptr;
shake_bond_t* d_shake_bonds = nullptr;
coord_t* d_coords = nullptr;
coord_t* d_xcoords = nullptr;
double* d_winv = nullptr;
int* d_total_iterations = nullptr;
int* d_mol_shake_offset = nullptr;
}  // namespace CudaShakeConstraints

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

int calc_shake_constraints_host() {
    using namespace CudaShakeConstraints;

    if (!is_initialized) {
        check_cudaMalloc((void**)&d_mol_n_shakes, sizeof(int) * n_molecules);
        check_cudaMalloc((void**)&d_shake_bonds, sizeof(shake_bond_t) * n_shake_constraints);
        check_cudaMalloc((void**)&d_coords, sizeof(coord_t) * n_atoms);
        check_cudaMalloc((void**)&d_xcoords, sizeof(coord_t) * n_atoms);
        check_cudaMalloc((void**)&d_winv, sizeof(double) * n_atoms);
        check_cudaMalloc((void**)&d_total_iterations, sizeof(int));
        check_cudaMalloc((void**)&d_mol_shake_offset, sizeof(int) * n_molecules);

        cudaMemcpy(d_mol_n_shakes, mol_n_shakes, sizeof(int) * n_molecules, cudaMemcpyHostToDevice);
        cudaMemcpy(d_shake_bonds, shake_bonds, sizeof(shake_bond_t) * n_shake_constraints, cudaMemcpyHostToDevice);
        cudaMemcpy(d_winv, winv, sizeof(double) * n_atoms, cudaMemcpyHostToDevice);

        int* mol_shake_offset_host = (int*)malloc(sizeof(int) * n_molecules);
        mol_shake_offset_host[0] = 0;
        for (int i = 1; i < n_molecules; i++) {
            mol_shake_offset_host[i] = mol_shake_offset_host[i - 1] + mol_n_shakes[i - 1];
        }
        cudaMemcpy(d_mol_shake_offset, mol_shake_offset_host, sizeof(int) * n_molecules, cudaMemcpyHostToDevice);

        free(mol_shake_offset_host);

        is_initialized = true;
    }

    cudaMemcpy(d_coords, coords, sizeof(coord_t) * n_atoms, cudaMemcpyHostToDevice);
    cudaMemcpy(d_xcoords, xcoords, sizeof(coord_t) * n_atoms, cudaMemcpyHostToDevice);
    int total_iterations_host = 0;
    cudaMemcpy(d_total_iterations, &total_iterations_host, sizeof(int), cudaMemcpyHostToDevice);

    int blocks = n_molecules;
    int threads = 32;
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

void cleanup_shake_constraints() {
    using namespace CudaShakeConstraints;

    if (is_initialized) {
        cudaFree(d_mol_n_shakes);
        cudaFree(d_shake_bonds);
        cudaFree(d_coords);
        cudaFree(d_xcoords);
        cudaFree(d_winv);
        cudaFree(d_total_iterations);
        cudaFree(d_mol_shake_offset);
        is_initialized = false;
    }
}