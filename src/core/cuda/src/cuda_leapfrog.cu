#include <iostream>

#include "common/include/context.h"
#include "cuda/include/cuda_force_accum.cuh"
#include "cuda/include/cuda_leapfrog.cuh"
#include "cuda/include/cuda_shake_constraints.cuh"
#include "cuda_utility.cuh"

namespace CudaLeapfrog {
}

__global__ void calc_leapfrog_kernel(
    atype_t* atypes,
    catype_t* catypes,
    vel_t* velocities,
    cuda_dvel_t* dvelocities,
    coord_t* coords,
    coord_t* xcoords,
    int n_atoms,
    int n_atoms_solute,
    real_t Tscale_solute,
    real_t Tscale_solvent,
    real_t dt) {
    int idx = blockIdx.x * blockDim.x + threadIdx.x;
    if (idx >= n_atoms) return;
    int i = idx;

    // Kernel implementation goes here
    real_t mass_i, winv_i;

    mass_i = catypes[atypes[i].code - 1].m;

    winv_i = 1 / mass_i;
    real_t scale = (i < n_atoms_solute) ? Tscale_solute : Tscale_solvent;
    const real_t dvel_x = force_accum_component_value(dvelocities[i].x);
    const real_t dvel_y = force_accum_component_value(dvelocities[i].y);
    const real_t dvel_z = force_accum_component_value(dvelocities[i].z);
    velocities[i].x = (velocities[i].x - dvel_x * dt * winv_i) * scale;
    velocities[i].y = (velocities[i].y - dvel_y * dt * winv_i) * scale;
    velocities[i].z = (velocities[i].z - dvel_z * dt * winv_i) * scale;

    xcoords[i].x = coords[i].x;
    xcoords[i].y = coords[i].y;
    xcoords[i].z = coords[i].z;

    coords[i].x += velocities[i].x * dt;
    coords[i].y += velocities[i].y * dt;
    coords[i].z += velocities[i].z * dt;
}

void calc_leapfrog_host() {
    auto& host = Context::instance();
    auto d_atypes = host.atypes->gpu_data_p;
    auto d_catypes = host.catypes->gpu_data_p;
    auto d_velocities = host.velocities->gpu_data_p;
    auto d_dvelocities = cuda_force_accum_buffer(host);
    auto d_coords = host.coords->gpu_data_p;
    auto d_xcoords = host.xcoords->gpu_data_p;

    int blockSize = 256;
    int numBlocks = (host.n_atoms + blockSize - 1) / blockSize;
    calc_leapfrog_kernel<<<numBlocks, blockSize>>>(
        d_atypes,
        d_catypes,
        d_velocities,
        d_dvelocities,
        d_coords,
        d_xcoords,
        host.n_atoms,
        host.n_atoms_solute,
        host.Tscale_solute,
        host.Tscale_solvent,
        host.dt);
    check_cuda(cudaDeviceSynchronize());

    host.velocities->download();
    host.coords->download();
    host.xcoords->download();

    // shake
    printf("n_shake_constraints: %d\n", host.n_shake_constraints);
    if (host.n_shake_constraints > 0) {
        calc_shake_constraints_host();
        auto& velocities = host.velocities->cpu_data_p;
        auto& coords = host.coords->cpu_data_p;
        auto* xcoords = host.xcoords->cpu_data_p;
        for (int i = 0; i < host.n_atoms; i++) {
            velocities[i].x = (coords[i].x - xcoords[i].x) / host.dt;
            velocities[i].y = (coords[i].y - xcoords[i].y) / host.dt;
            velocities[i].z = (coords[i].z - xcoords[i].z) / host.dt;
        }
        host.velocities->upload();
    }
}

void init_leapfrog_kernel_data() {}

void cleanup_leapfrog() {}
