#include <iostream>

#include "cuda/include/CudaContext.cuh"
#include "cuda/include/CudaLeapfrog.cuh"
#include "cuda/include/CudaShakeConstraints.cuh"
#include "utils.h"

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
    CudaContext& ctx = CudaContext::instance();
    ctx.sync_all_to_device();
    auto d_atypes = ctx.d_atypes;
    auto d_catypes = ctx.d_catypes;
    auto d_velocities = ctx.d_velocities;
    auto d_dvelocities = ctx.d_dvelocities;
    auto d_coords = ctx.d_coords;
    auto d_xcoords = ctx.d_xcoords;

    int blockSize = 256;
    int numBlocks = (n_atoms + blockSize - 1) / blockSize;
    calc_leapfrog_kernel<<<numBlocks, blockSize>>>(
        d_atypes,
        d_catypes,
        d_velocities,
        d_dvelocities,
        d_coords,
        d_xcoords,
        n_atoms,
        n_atoms_solute,
        Tscale_solute,
        Tscale_solvent,
        dt);
    cudaDeviceSynchronize();

    cudaMemcpy(velocities, d_velocities, sizeof(vel_t) * n_atoms, cudaMemcpyDeviceToHost);
    cudaMemcpy(dvelocities, d_dvelocities, sizeof(dvel_t) * n_atoms, cudaMemcpyDeviceToHost);
    cudaMemcpy(coords, d_coords, sizeof(coord_t) * n_atoms, cudaMemcpyDeviceToHost);
    cudaMemcpy(xcoords, d_xcoords, sizeof(coord_t) * n_atoms, cudaMemcpyDeviceToHost);

    // shake
    printf("n_shake_constraints: %d\n", n_shake_constraints);
    if (n_shake_constraints > 0) {
        calc_shake_constraints_host();
        for (int i = 0; i < n_atoms; i++) {
            velocities[i].x = (coords[i].x - xcoords[i].x) / dt;
            velocities[i].y = (coords[i].y - xcoords[i].y) / dt;
            velocities[i].z = (coords[i].z - xcoords[i].z) / dt;
        }
    }
}