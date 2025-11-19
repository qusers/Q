#include "cuda/include/CudaBondForce.cuh"
#include "utils.h"

namespace CudaBondForce {
bool is_initialized = false;
bond_t* d_bonds = nullptr;
coord_t* d_coords = nullptr;
cbond_t* d_cbonds = nullptr;
dvel_t* d_dvelocities = nullptr;
double* d_energy_sum = nullptr;
}  // namespace CudaBondForce

__global__ void calc_bond_forces_kernel(int start, int end, bond_t* bonds, coord_t* coords, cbond_t* cbonds, dvel_t* dvelocities, double* energy_sum) {
    int idx = blockIdx.x * blockDim.x + threadIdx.x + start;
    if (idx >= end) return;
    bond_t bond = bonds[idx];
    coord_t ri = coords[bond.ai - 1];
    coord_t rj = coords[bond.aj - 1];
    cbond_t cbond = cbonds[bond.code - 1];

    double dx = rj.x - ri.x;
    double dy = rj.y - ri.y;
    double dz = rj.z - ri.z;
    double r = sqrt(dx * dx + dy * dy + dz * dz);

    double dr = r - cbond.b0;
    double energy = 0.5 * cbond.kb * dr * dr;

    atomicAdd(energy_sum, energy);

    // update forces
    double f = cbond.kb * dr / r;
    atomicAdd(&dvelocities[bond.aj - 1].x, f * dx);
    atomicAdd(&dvelocities[bond.aj - 1].y, f * dy);
    atomicAdd(&dvelocities[bond.aj - 1].z, f * dz);
    atomicAdd(&dvelocities[bond.ai - 1].x, -f * dx);
    atomicAdd(&dvelocities[bond.ai - 1].y, -f * dy);
    atomicAdd(&dvelocities[bond.ai - 1].z, -f * dz);
}

double calc_bond_forces_host(int start, int end) {
    using namespace CudaBondForce;
    int N = end - start;
    if (N <= 0) return 0.0;
    int blockSize = 256;
    int numBlocks = (N + blockSize - 1) / blockSize;

    if (!is_initialized) {
        check_cudaMalloc((void**)&d_bonds, sizeof(bond_t) * n_bonds);
        check_cudaMalloc((void**)&d_coords, sizeof(coord_t) * n_atoms);
        check_cudaMalloc((void**)&d_cbonds, sizeof(cbond_t) * n_cbonds);
        check_cudaMalloc((void**)&d_dvelocities, sizeof(dvel_t) * n_atoms);
        check_cudaMalloc((void**)&d_energy_sum, sizeof(double));
        is_initialized = true;
    }

    cudaMemcpy(d_bonds, bonds, sizeof(bond_t) * n_bonds, cudaMemcpyHostToDevice);
    cudaMemcpy(d_coords, coords, sizeof(coord_t) * n_atoms, cudaMemcpyHostToDevice);
    cudaMemcpy(d_cbonds, cbonds, sizeof(cbond_t) * n_cbonds, cudaMemcpyHostToDevice);
    cudaMemcpy(d_dvelocities, dvelocities, sizeof(dvel_t) * n_atoms, cudaMemcpyHostToDevice);
    double zero = 0.0;
    cudaMemcpy(d_energy_sum, &zero, sizeof(double), cudaMemcpyHostToDevice);
    calc_bond_forces_kernel<<<numBlocks, blockSize>>>(start, end, d_bonds, d_coords, d_cbonds, d_dvelocities, d_energy_sum);
    cudaDeviceSynchronize();
    cudaMemcpy(&zero, d_energy_sum, sizeof(double), cudaMemcpyDeviceToHost);
    cudaMemcpy(dvelocities, d_dvelocities, sizeof(dvel_t) * n_atoms, cudaMemcpyDeviceToHost);
    return zero;
}

void cleanup_bond_force() {
    using namespace CudaBondForce;
    if (is_initialized) {
        cudaFree(d_bonds);
        cudaFree(d_coords);
        cudaFree(d_cbonds);
        cudaFree(d_dvelocities);
        cudaFree(d_energy_sum);
        is_initialized = false;
    }
}