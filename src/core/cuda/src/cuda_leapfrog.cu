#include <iostream>

#include "common/include/context.h"
#include "cuda/include/cuda_context.cuh"
#include "cuda/include/cuda_leapfrog.cuh"
#include "cuda/include/cuda_shake_constraints.cuh"
#include "cuda_utility.cuh"

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
    double Tscale_solute,
    double Tscale_solvent,
    double dt) {
    int idx = blockIdx.x * blockDim.x + threadIdx.x;
    if (idx >= n_atoms) return;
    int i = idx;

    // Kernel implementation goes here
    double mass_i, winv_i;

    mass_i = catypes[atypes[i].code - 1].m;

    winv_i = 1 / mass_i;
    double scale = (i < n_atoms_solute) ? Tscale_solute : Tscale_solvent;
    velocities[i].x = (velocities[i].x - dvelocities[i].x * dt * winv_i) * scale;
    velocities[i].y = (velocities[i].y - dvelocities[i].y * dt * winv_i) * scale;
    velocities[i].z = (velocities[i].z - dvelocities[i].z * dt * winv_i) * scale;

    xcoords[i].x = coords[i].x;
    xcoords[i].y = coords[i].y;
    xcoords[i].z = coords[i].z;

    coords[i].x += velocities[i].x * dt;
    coords[i].y += velocities[i].y * dt;
    coords[i].z += velocities[i].z * dt;
}

void calc_leapfrog_host() {
    auto& host = Context::instance();
    auto& ctx = CudaContext::instance();
    auto d_atypes = host.atypes->gpu_data_p;
    auto d_catypes = host.catypes->gpu_data_p;
    auto d_velocities = ctx.d_velocities;
    auto d_dvelocities = ctx.d_dvelocities;
    auto d_coords = ctx.d_coords;
    auto d_xcoords = ctx.d_xcoords;

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

    check_cuda(cudaMemcpy(host.velocities.data(), d_velocities, sizeof(vel_t) * host.n_atoms, cudaMemcpyDeviceToHost));
    check_cuda(cudaMemcpy(host.dvelocities.data(), d_dvelocities, sizeof(dvel_t) * host.n_atoms, cudaMemcpyDeviceToHost));
    check_cuda(cudaMemcpy(host.coords.data(), d_coords, sizeof(coord_t) * host.n_atoms, cudaMemcpyDeviceToHost));
    check_cuda(cudaMemcpy(host.xcoords.data(), d_xcoords, sizeof(coord_t) * host.n_atoms, cudaMemcpyDeviceToHost));

    // shake
    // todo: Here is some problem, it writes into cpu memory, but we use gpu..
    printf("n_shake_constraints: %d\n", host.n_shake_constraints);
    if (host.n_shake_constraints > 0) {
        calc_shake_constraints_host();
        for (int i = 0; i < host.n_atoms; i++) {
            host.velocities[i].x = (host.coords[i].x - host.xcoords[i].x) / host.dt;
            host.velocities[i].y = (host.coords[i].y - host.xcoords[i].y) / host.dt;
            host.velocities[i].z = (host.coords[i].z - host.xcoords[i].z) / host.dt;
        }
    }
}

void init_leapfrog_kernel_data() {}

void cleanup_leapfrog() {}
