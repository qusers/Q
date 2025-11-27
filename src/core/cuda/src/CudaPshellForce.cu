#include "cuda/include/CudaContext.cuh"
#include "cuda/include/CudaPshellForce.cuh"
#include "utils.h"

__global__ void calc_pshell_force_kernel(
    int n_atoms_solute,
    bool* shell,
    bool* excluded,
    coord_t* coords,
    coord_t* coords_top,
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
        dr.x = coords[i].x - coords_top[i].x;
        dr.y = coords[i].y - coords_top[i].y;
        dr.z = coords[i].z - coords_top[i].z;
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
    CudaContext& ctx = CudaContext::instance();
    ctx.sync_all_to_device();

    auto d_shell = ctx.d_shell;
    auto d_excluded = ctx.d_excluded;
    auto d_coords = ctx.d_coords;
    auto d_coords_top = ctx.d_coords_top;
    auto d_dvelocities = ctx.d_dvelocities;

    double* d_ufix_energy;
    double* d_ushell_energy;
    check_cudaMalloc((void**)&d_ufix_energy, sizeof(double));
    check_cudaMalloc((void**)&d_ushell_energy, sizeof(double));
    cudaMemset(d_ufix_energy, 0, sizeof(double));
    cudaMemset(d_ushell_energy, 0, sizeof(double));

    int blockSize = 256;
    int numBlocks = (n_atoms_solute + blockSize - 1) / blockSize;
    calc_pshell_force_kernel<<<numBlocks, blockSize>>>(
        n_atoms_solute,
        d_shell,
        d_excluded,
        d_coords,
        d_coords_top,
        d_ufix_energy,
        d_ushell_energy,
        d_dvelocities);
    cudaDeviceSynchronize();
    double ufix_energy;
    double ushell_energy;
    cudaMemcpy(&ufix_energy, d_ufix_energy, sizeof(double), cudaMemcpyDeviceToHost);
    cudaMemcpy(&ushell_energy, d_ushell_energy, sizeof(double), cudaMemcpyDeviceToHost);
    cudaMemcpy(dvelocities, d_dvelocities, sizeof(dvel_t) * n_atoms, cudaMemcpyDeviceToHost);

    E_restraint.Ufix += ufix_energy;
    E_restraint.Ushell += ushell_energy;
    // ctx.sync_all_to_host();
}