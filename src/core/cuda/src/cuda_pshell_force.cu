#include "cuda/include/cuda_pshell_force.cuh"
#include "common/include/constants.h"
#include "common/include/context.h"
#include "cuda_utility.cuh"
#include <iostream>
namespace CudaPshellForce {
bool is_initialized = false;
double* d_ufix_energy;
double* d_ushell_energy;

}  // namespace CudaPshellForce
__global__ void calc_pshell_force_kernel(
    int n_atoms_solute,
    bool* shell,
    bool* excluded,
    coord_t* coords,
    coord_t* coords_init,
    double* ufix_energy,
    double* ushell_energy,
    dvel_t* dvelocities) {
    int i = blockIdx.x * blockDim.x + threadIdx.x;
    if (i >= n_atoms_solute) return;

    coord_t dr;
    double k, r2, ener;

    if (shell[i] || excluded[i]) {
        // printf("i = %d excluded = %s shell = %s\n", i, excluded[i] ? "True" : "False", shell[i] ? "True" : "False");
        if (excluded[i]) {
            k = k_fix;
        } else {
            k = k_pshell;
        }
        dr.x = coords[i].x - coords_init[i].x;
        dr.y = coords[i].y - coords_init[i].y;
        dr.z = coords[i].z - coords_init[i].z;
        r2 = pow(dr.x, 2) + pow(dr.y, 2) + pow(dr.z, 2);
        ener = 0.5 * k * r2;
        // printf("dr = %f %f %f\n", dr.x, dr.y, dr.z);

        if (excluded[i]) atomicAdd(ufix_energy, ener);
        if (shell[i]) atomicAdd(ushell_energy, ener);

        atomicAdd(&dvelocities[i].x, k * dr.x);
        atomicAdd(&dvelocities[i].y, k * dr.y);
        atomicAdd(&dvelocities[i].z, k * dr.z);
    }
}

void calc_pshell_forces_host() {
    auto& host = Context::instance();
    using namespace CudaPshellForce;

    auto d_shell = host.shell->gpu_data_p;
    auto d_excluded = host.excluded->gpu_data_p;
    auto d_coords = host.coords->gpu_data_p;
    auto d_coords_init = host.coords_init->gpu_data_p;
    auto d_dvelocities = host.dvelocities->gpu_data_p;

    cudaMemset(d_ufix_energy, 0, sizeof(double));
    cudaMemset(d_ushell_energy, 0, sizeof(double));

    int blockSize = 256;
    int numBlocks = (host.n_atoms_solute + blockSize - 1) / blockSize;
    calc_pshell_force_kernel<<<numBlocks, blockSize>>>(
        host.n_atoms_solute,
        d_shell,
        d_excluded,
        d_coords,
        d_coords_init,
        d_ufix_energy,
        d_ushell_energy,
        d_dvelocities);
    cudaDeviceSynchronize();
    double ufix_energy;
    double ushell_energy;
    cudaMemcpy(&ufix_energy, d_ufix_energy, sizeof(double), cudaMemcpyDeviceToHost);
    cudaMemcpy(&ushell_energy, d_ushell_energy, sizeof(double), cudaMemcpyDeviceToHost);

    host.E_restraint.Ufix += ufix_energy;
    host.E_restraint.Ushell += ushell_energy;
    // ctx.sync_all_to_host();
}

void init_pshell_force_kernel_data() {
    using namespace CudaPshellForce;
    if (!is_initialized) {
        check_cudaMalloc((void**)&d_ufix_energy, sizeof(double));
        check_cudaMalloc((void**)&d_ushell_energy, sizeof(double));
        is_initialized = true;
    }
}

void cleanup_pshell_force() {
    using namespace CudaPshellForce;
    if (is_initialized) {
        cudaFree(d_ufix_energy);
        cudaFree(d_ushell_energy);
        is_initialized = false;
    }
}
