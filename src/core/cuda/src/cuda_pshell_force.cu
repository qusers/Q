#include <iostream>

#include "common/include/constants.h"
#include "common/include/context.h"
#include "cuda/include/cuda_pshell_force.cuh"
#include "cuda_force_accumulation.cuh"
#include "cuda_utility.cuh"
namespace CudaPshellForce {
}  // namespace CudaPshellForce
__global__ void calc_pshell_force_kernel(
    int n_atoms_solute,
    bool* shell,
    bool* excluded,
    coord_t* coords,
    coord_t* coords_init,
    energy_accum_t* energy,
    dvel_t* dvelocities) {
    int i = blockIdx.x * blockDim.x + threadIdx.x;
    if (i >= n_atoms_solute) return;

    if (shell[i] || excluded[i]) {
        // printf("i = %d excluded = %s shell = %s\n", i, excluded[i] ? "True" : "False", shell[i] ? "True" : "False");
        const double k = excluded[i] ? k_fix : k_pshell;
        const double dx = coords[i].x - coords_init[i].x;
        const double dy = coords[i].y - coords_init[i].y;
        const double dz = coords[i].z - coords_init[i].z;
        const double r2 = dx * dx + dy * dy + dz * dz;
        const double ener = 0.5 * k * r2;
        // printf("dr = %f %f %f\n", dr.x, dr.y, dr.z);

        if (excluded[i]) atomic_add_energy(&energy[E_RESTR_FIX], ener);
        if (shell[i]) atomic_add_energy(&energy[E_RESTR_SHELL], ener);

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

    int blockSize = 256;
    int numBlocks = (host.n_atoms_solute + blockSize - 1) / blockSize;
    calc_pshell_force_kernel<<<numBlocks, blockSize>>>(
        host.n_atoms_solute,
        d_shell,
        d_excluded,
        d_coords,
        d_coords_init,
        host.energy.device(),
        d_dvelocities);
}

void init_pshell_force_kernel_data() {
}

void cleanup_pshell_force() {
}
