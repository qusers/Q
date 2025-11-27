#include "cuda/include/CudaBondForce.cuh"
#include "cuda/include/CudaContext.cuh"
#include "utils.h"

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
    int N = end - start;
    if (N <= 0) return 0.0;
    int blockSize = 256;
    int numBlocks = (N + blockSize - 1) / blockSize;

    double energy = 0.0;
    double* d_energy_sum;
    check_cudaMalloc((void**)&d_energy_sum, sizeof(double));
    cudaMemcpy(d_energy_sum, &energy, sizeof(double), cudaMemcpyHostToDevice);

    CudaContext& context = CudaContext::instance();
    context.sync_all_to_device();
    bond_t* d_bonds = context.d_bonds;
    coord_t* d_coords = context.d_coords;
    cbond_t* d_cbonds = context.d_cbonds;
    dvel_t* d_dvelocities = context.d_dvelocities;

    calc_bond_forces_kernel<<<numBlocks, blockSize>>>(start, end, d_bonds, d_coords, d_cbonds, d_dvelocities, d_energy_sum);
    cudaDeviceSynchronize();
    cudaMemcpy(&energy, d_energy_sum, sizeof(double), cudaMemcpyDeviceToHost);
    cudaMemcpy(dvelocities, d_dvelocities, sizeof(dvel_t) * n_atoms, cudaMemcpyDeviceToHost);

    cudaFree(d_energy_sum);
    return energy;
}
