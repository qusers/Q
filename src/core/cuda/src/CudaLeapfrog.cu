#include "cuda/include/CudaLeapfrog.cuh"
#include "utils.h"

namespace CUDALeapfrog {
bool is_initialized = false;
atype_t* d_atypes = nullptr;
catype_t* d_catypes = nullptr;
vel_t* d_velocities = nullptr;
dvel_t* d_dvelocities = nullptr;
coord_t* d_coords = nullptr;
coord_t* d_xcoords = nullptr;

}  // namespace CUDALeapfrog





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
    double dt
) {
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
    using namespace CUDALeapfrog;
    if (!is_initialized) {
        check_cudaMalloc((void**)&d_atypes, sizeof(atype_t) * n_atypes);
        check_cudaMalloc((void**)&d_catypes, sizeof(catype_t) * n_catypes);
        check_cudaMalloc((void**)&d_velocities, sizeof(vel_t) * n_atoms);
        check_cudaMalloc((void**)&d_dvelocities, sizeof(dvel_t) * n_atoms);
        check_cudaMalloc((void**)&d_coords, sizeof(coord_t) * n_atoms);
        check_cudaMalloc((void**)&d_xcoords, sizeof(coord_t) * n_atoms);

        cudaMemcpy(d_atypes, atypes, sizeof(atype_t) * n_atypes, cudaMemcpyHostToDevice);
        cudaMemcpy(d_catypes, catypes, sizeof(catype_t) * n_catypes, cudaMemcpyHostToDevice);

        is_initialized = true;
    }
    cudaMemcpy(d_velocities, velocities, sizeof(vel_t) * n_atoms, cudaMemcpyHostToDevice);
    cudaMemcpy(d_dvelocities, dvelocities, sizeof(dvel_t) * n_atoms, cudaMemcpyHostToDevice);
    cudaMemcpy(d_coords, coords, sizeof(coord_t) * n_atoms, cudaMemcpyHostToDevice);

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
        dt
    );
    cudaDeviceSynchronize();
    cudaMemcpy(velocities, d_velocities, sizeof(vel_t) * n_atoms, cudaMemcpyDeviceToHost);
    cudaMemcpy(dvelocities, d_dvelocities, sizeof(dvel_t) * n_atoms, cudaMemcpyDeviceToHost);
    cudaMemcpy(coords, d_coords, sizeof(coord_t) * n_atoms, cudaMemcpyDeviceToHost);
    cudaMemcpy(xcoords, d_xcoords, sizeof(coord_t) * n_atoms, cudaMemcpyDeviceToHost);

    
}
void cleanup_leapfrog() {
}
