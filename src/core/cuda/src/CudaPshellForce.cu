#include "cuda/include/CudaPshellForce.cuh"
#include "utils.h"
namespace CudaPshellForce {
bool is_initialized = false;
bool* d_shell;
bool* d_excluded;
coord_t* d_coords;
coord_t* d_coords_top;
double* d_ufix_energy;
double* d_ushell_energy;
dvel_t* d_dvelocities;
}  // namespace CudaPshellForce

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
    using namespace CudaPshellForce;
    if (!is_initialized) {
        check_cudaMalloc((void**)&d_shell, sizeof(bool) * n_atoms);
        check_cudaMalloc((void**)&d_excluded, sizeof(bool) * n_atoms);
        check_cudaMalloc((void**)&d_coords, sizeof(coord_t) * n_atoms);
        check_cudaMalloc((void**)&d_coords_top, sizeof(coord_t) * n_atoms);
        check_cudaMalloc((void**)&d_ufix_energy, sizeof(double));
        check_cudaMalloc((void**)&d_ushell_energy, sizeof(double));
        check_cudaMalloc((void**)&d_dvelocities, sizeof(dvel_t) * n_atoms);

        cudaMemcpy(d_shell, shell, sizeof(bool) * n_atoms, cudaMemcpyHostToDevice);
        cudaMemcpy(d_excluded, excluded, sizeof(bool) * n_atoms, cudaMemcpyHostToDevice);
        is_initialized = true;
    }

    cudaMemcpy(d_coords, coords, sizeof(coord_t) * n_atoms, cudaMemcpyHostToDevice);
    cudaMemcpy(d_coords_top, coords_top, sizeof(coord_t) * n_atoms, cudaMemcpyHostToDevice);
    cudaMemcpy(d_dvelocities, dvelocities, sizeof(dvel_t) * n_atoms, cudaMemcpyHostToDevice);

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
}

void cleanup_pshell_force() {
    using namespace CudaPshellForce;
    if (is_initialized) {
        cudaFree(d_shell);
        cudaFree(d_excluded);
        cudaFree(d_coords);
        cudaFree(d_coords_top);
        cudaFree(d_ufix_energy);
        cudaFree(d_ushell_energy);
        cudaFree(d_dvelocities);
        is_initialized = false;
    }

}