#include "cuda/include/cuda_bond_force.cuh"
#include "context.h"
#include "cuda_utility.cuh"
namespace CudaBondForce {
bool is_initialized = false;
double* d_energy_sum;
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
    int N = end - start;
    if (N <= 0) return 0.0;
    using namespace CudaBondForce;
    int blockSize = 256;
    int numBlocks = (N + blockSize - 1) / blockSize;

    double energy = 0.0;
    cudaMemcpy(d_energy_sum, &energy, sizeof(double), cudaMemcpyHostToDevice);

    auto& host_ctx = Context::instance();
    bond_t* d_bonds = host_ctx.bonds->gpu_data_p;
    coord_t* d_coords = host_ctx.coords->gpu_data_p;
    cbond_t* d_cbonds = host_ctx.cbonds->gpu_data_p;
    dvel_t* d_dvelocities = host_ctx.dvelocities->gpu_data_p;

    calc_bond_forces_kernel<<<numBlocks, blockSize>>>(start, end, d_bonds, d_coords, d_cbonds, d_dvelocities, d_energy_sum);
    cudaDeviceSynchronize();
    cudaMemcpy(&energy, d_energy_sum, sizeof(double), cudaMemcpyDeviceToHost);

    return energy;
}

void init_bond_force_kernel_data() {
    using namespace CudaBondForce;
    if (!is_initialized) {
        check_cudaMalloc((void**)&d_energy_sum, sizeof(double));
        is_initialized = true;
    }
}

void cleanup_bond_force() {
    using namespace CudaBondForce;
    if (is_initialized) {
        cudaFree(d_energy_sum);
        is_initialized = false;
    }
}
