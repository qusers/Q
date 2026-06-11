#include "cuda/include/cuda_pshell_force.cuh"
#include "common/include/constants.h"
#include "common/include/context.h"
#include "cuda_utility.cuh"
#include <iostream>
#include "cuda_force_accumulation.cuh"
namespace CudaPshellForce {
bool is_initialized = false;
energy_accum_t* d_ufix_energy;
energy_accum_t* d_ushell_energy;

}  // namespace CudaPshellForce
__global__ void calc_pshell_force_kernel(
    int n_atoms_solute,
    bool* shell,
    bool* excluded,
    coord_t* coords,
    coord_t* coords_init,
    energy_accum_t* ufix_energy,
    energy_accum_t* ushell_energy,
    dvel_t* dvelocities) {
    int i = blockIdx.x * blockDim.x + threadIdx.x;
    if (i >= n_atoms_solute) return;

    if (shell[i] || excluded[i]) {
        // printf("i = %d excluded = %s shell = %s\n", i, excluded[i] ? "True" : "False", shell[i] ? "True" : "False");
        const double k = excluded[i] ? k_fix : k_pshell;
        const double dx = static_cast<double>(coords[i].x) - static_cast<double>(coords_init[i].x);
        const double dy = static_cast<double>(coords[i].y) - static_cast<double>(coords_init[i].y);
        const double dz = static_cast<double>(coords[i].z) - static_cast<double>(coords_init[i].z);
        const double r2 = dx * dx + dy * dy + dz * dz;
        const double ener = 0.5 * k * r2;
        // printf("dr = %f %f %f\n", dr.x, dr.y, dr.z);

        if (excluded[i]) atomic_add_energy(ufix_energy, ener);
        if (shell[i]) atomic_add_energy(ushell_energy, ener);

        atomic_add_force(&dvelocities[i].x, k * dx);
        atomic_add_force(&dvelocities[i].y, k * dy);
        atomic_add_force(&dvelocities[i].z, k * dz);
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

    cudaMemset(d_ufix_energy, 0, sizeof(energy_accum_t));
    cudaMemset(d_ushell_energy, 0, sizeof(energy_accum_t));

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
    energy_accum_t ufix_energy;
    energy_accum_t ushell_energy;
    cudaMemcpy(&ufix_energy, d_ufix_energy, sizeof(energy_accum_t), cudaMemcpyDeviceToHost);
    cudaMemcpy(&ushell_energy, d_ushell_energy, sizeof(energy_accum_t), cudaMemcpyDeviceToHost);

    host.E_restraint.Ufix += energy_from_accum(ufix_energy);
    host.E_restraint.Ushell += energy_from_accum(ushell_energy);
    // ctx.sync_all_to_host();
}

void init_pshell_force_kernel_data() {
    using namespace CudaPshellForce;
    if (!is_initialized) {
        check_cudaMalloc((void**)&d_ufix_energy, sizeof(energy_accum_t));
        check_cudaMalloc((void**)&d_ushell_energy, sizeof(energy_accum_t));
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
