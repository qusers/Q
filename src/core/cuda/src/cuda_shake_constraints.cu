#include <cstdio>
#include <cstdlib>
#include <iostream>

#include "cuda/include/cuda_shake_constraints.cuh"
#include "common/include/constants.h"
#include "common/include/context.h"

namespace CudaShakeConstraints {

bool is_initialized = false;
int* d_mol_shake_offset;
int* d_total_iterations;
int* d_failed;
}  // namespace CudaShakeConstraints

__global__ void calc_shake_constraints_kernel(
    int n_molecules,
    int* mol_n_shakes,
    int* mol_n_shake_groups,
    int* mol_shake_group_offset,
    int* shake_group_offset,
    int* shake_group_size,
    int* shake_group_indices,
    shake_bond_t* shake_bonds,
    coord_t* coords,
    coord_t* xcoords,
    real_t* winv,
    int* total_iterations,
    int* failed,
    int* mol_shake_offset) {
    int idx = blockIdx.x;
    if (idx >= n_molecules) return;

    int mol = idx;
    int mol_shake_count = mol_n_shakes[mol];
    if (mol_shake_count == 0) return;

    int ai, aj;
    real_t xij2, diff, corr, scp, xxij2;
    coord_t xij, xxij;
    __shared__ int mol_converged;

    int shake = mol_shake_offset[mol];
    int group_base = mol_shake_group_offset[mol];
    int group_count = mol_n_shake_groups[mol];
    int n_iterations = 0;

    for (int i = threadIdx.x; i < mol_shake_count; i += blockDim.x) {
        shake_bonds[shake + i].ready = false;
    }
    __syncthreads();

    do {
        for (int group_i = 0; group_i < group_count; group_i++) {
            int group = group_base + group_i;
            int group_offset = shake_group_offset[group];
            int group_size = shake_group_size[group];

            for (int local_i = threadIdx.x; local_i < group_size; local_i += blockDim.x) {
                int shake_i = shake_group_indices[group_offset + local_i];
                if (!shake_bonds[shake_i].ready) {
                    ai = shake_bonds[shake_i].ai - 1;
                    aj = shake_bonds[shake_i].aj - 1;

                    xij.x = coords[ai].x - coords[aj].x;
                    xij.y = coords[ai].y - coords[aj].y;
                    xij.z = coords[ai].z - coords[aj].z;
                    xij2 = xij.x * xij.x + xij.y * xij.y + xij.z * xij.z;
                    diff = shake_bonds[shake_i].dist2 - xij2;
                    if (fabs(diff) < shake_tol * shake_bonds[shake_i].dist2) {
                        shake_bonds[shake_i].ready = true;
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
            __syncthreads();
        }

        n_iterations++;

        if (threadIdx.x == 0) {
            mol_converged = 1;
            for (int i = 0; i < mol_shake_count; i++) {
                if (!shake_bonds[shake + i].ready) {
                    mol_converged = 0;
                    break;
                }
            }
        }
        __syncthreads();
    } while (n_iterations < shake_max_iter && mol_converged == 0);

    if (threadIdx.x == 0) {
        if (mol_converged == 0) {
            for (int i = 0; i < mol_shake_count; i++) {
                if (shake_bonds[shake + i].ready) continue;
                ai = shake_bonds[shake + i].ai - 1;
                aj = shake_bonds[shake + i].aj - 1;

                xxij.x = xcoords[ai].x - xcoords[aj].x;
                xxij.y = xcoords[ai].y - xcoords[aj].y;
                xxij.z = xcoords[ai].z - xcoords[aj].z;
                xxij2 = xxij.x * xxij.x + xxij.y * xxij.y + xxij.z * xxij.z;
                printf(">>> Shake failed, i = %d,j = %d, d = %f, d0 = %f\n", ai + 1, aj + 1, sqrt(xxij2), sqrt(shake_bonds[shake + i].dist2));
            }
            atomicExch(failed, 1);
            return;
        }

        atomicAdd(total_iterations, n_iterations);
    }
}

void init_shake_constraints_kernel_data() {
    auto& host = Context::instance();
    auto *mol_n_shakes = host.mol_n_shakes->cpu_data_p;
    using namespace CudaShakeConstraints;
    if (!is_initialized) {
        if (host.n_molecules > 0) {
            int* mol_shake_offset_host = (int*)malloc(sizeof(int) * host.n_molecules);
            mol_shake_offset_host[0] = 0;
            for (int i = 1; i < host.n_molecules; i++) {
                mol_shake_offset_host[i] = mol_shake_offset_host[i - 1] + mol_n_shakes[i - 1];
            }
            check_cudaMalloc((void**)&d_mol_shake_offset, sizeof(int) * host.n_molecules);
            cudaMemcpy(d_mol_shake_offset, mol_shake_offset_host, sizeof(int) * host.n_molecules, cudaMemcpyHostToDevice);
            free(mol_shake_offset_host);
        }

        check_cudaMalloc((void**)&d_total_iterations, sizeof(int));
        check_cudaMalloc((void**)&d_failed, sizeof(int));

        is_initialized = true;
    }
}

void cleanup_shake_constraints() {
    using namespace CudaShakeConstraints;
    if (is_initialized) {
        cudaFree(d_mol_shake_offset);
        cudaFree(d_total_iterations);
        cudaFree(d_failed);
        is_initialized = false;
    }
}

int calc_shake_constraints_host() {
    auto& host = Context::instance();
    if (host.n_molecules == 0 || host.n_shake_constraints == 0) return 0;
    using namespace CudaShakeConstraints;
    int total_iterations_host = 0;
    int failed_host = 0;
    check_cuda(cudaMemcpy(d_total_iterations, &total_iterations_host, sizeof(int), cudaMemcpyHostToDevice));
    check_cuda(cudaMemcpy(d_failed, &failed_host, sizeof(int), cudaMemcpyHostToDevice));

    int blocks = host.n_molecules;
    int threads = 32;

    auto d_mol_n_shakes = host.mol_n_shakes->gpu_data_p;
    auto d_mol_n_shake_groups = host.mol_n_shake_groups->gpu_data_p;
    auto d_mol_shake_group_offset = host.mol_shake_group_offset->gpu_data_p;
    auto d_shake_group_offset = host.shake_group_offset->gpu_data_p;
    auto d_shake_group_size = host.shake_group_size->gpu_data_p;
    auto d_shake_group_indices = host.shake_group_indices->gpu_data_p;
    auto d_shake_bonds = host.shake_bonds->gpu_data_p;
    auto d_coords = host.coords->gpu_data_p;
    auto d_xcoords = host.xcoords->gpu_data_p;
    auto d_winv = host.winv->gpu_data_p;

    calc_shake_constraints_kernel<<<blocks, threads>>>(
        host.n_molecules,
        d_mol_n_shakes,
        d_mol_n_shake_groups,
        d_mol_shake_group_offset,
        d_shake_group_offset,
        d_shake_group_size,
        d_shake_group_indices,
        d_shake_bonds,
        d_coords,
        d_xcoords,
        d_winv,
        d_total_iterations,
        d_failed,
        d_mol_shake_offset);
    check_cuda(cudaDeviceSynchronize());
    check_cuda(cudaMemcpy(&failed_host, d_failed, sizeof(int), cudaMemcpyDeviceToHost));
    check_cuda(cudaMemcpy(&total_iterations_host, d_total_iterations, sizeof(int), cudaMemcpyDeviceToHost));
    host.coords->download();
    if (failed_host != 0) {
        std::printf(">>> FATAL: shake failure\n");
        std::exit(EXIT_FAILURE);
    }
    return host.n_molecules == 0 ? 0 : total_iterations_host / host.n_molecules;
}
