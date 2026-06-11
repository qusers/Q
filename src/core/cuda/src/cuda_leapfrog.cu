#include <iostream>

#include "common/include/context.h"
#include "cuda/include/cuda_leapfrog.cuh"
#include "cuda_utility.cuh"
#include "cuda_force_accumulation.cuh"

namespace CudaLeapfrog {

}

__global__ void calc_leapfrog_kernel(
    atype_t* atypes,
    catype_t* catypes,
    vel_t* velocities,
    dvel_t* dvelocities,
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

    const real_t fx = force_from_accum(dvelocities[i].x);
    const real_t fy = force_from_accum(dvelocities[i].y);
    const real_t fz = force_from_accum(dvelocities[i].z);

    velocities[i].x = (velocities[i].x - fx * dt * winv_i) * scale;
    velocities[i].y = (velocities[i].y - fy * dt * winv_i) * scale;
    velocities[i].z = (velocities[i].z - fz * dt * winv_i) * scale;

    xcoords[i].x = coords[i].x;
    xcoords[i].y = coords[i].y;
    xcoords[i].z = coords[i].z;

    coords[i].x += velocities[i].x * dt;
    coords[i].y += velocities[i].y * dt;
    coords[i].z += velocities[i].z * dt;
}

__global__ void update_velocities_from_positions_kernel(
    vel_t* velocities,
    const coord_t* coords,
    const coord_t* xcoords,
    int n_atoms,
    real_t dt) {
    int idx = blockIdx.x * blockDim.x + threadIdx.x;
    if (idx >= n_atoms) return;

    velocities[idx].x = (coords[idx].x - xcoords[idx].x) / dt;
    velocities[idx].y = (coords[idx].y - xcoords[idx].y) / dt;
    velocities[idx].z = (coords[idx].z - xcoords[idx].z) / dt;
}

void calc_leapfrog_host(Shake &shake) {
    auto& host = Context::instance();
    auto d_atypes = host.atypes->gpu_data_p;
    auto d_catypes = host.catypes->gpu_data_p;
    auto d_velocities = host.velocities->gpu_data_p;
    auto d_dvelocities = host.dvelocities->gpu_data_p;
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

    // shake
    if (shake.data().n_constraints > 0) {
        shake.apply(host);
        update_velocities_from_positions_kernel<<<numBlocks, blockSize>>>(
            d_velocities,
            d_coords,
            d_xcoords,
            host.n_atoms,
            host.dt);
    }
    check_cuda(cudaDeviceSynchronize());
}

void init_leapfrog_kernel_data() {}

void cleanup_leapfrog() {}
