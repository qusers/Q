#include "cuda/include/cuda_bond_force.cuh"
#include "context.h"
#include "cuda_utility.cuh"
#include "cuda_force_accumulation.cuh"
namespace CudaBondForce {
bool is_initialized = false;
energy_accum_t* d_energy_sum;
}  // namespace CudaBondForce
__global__ void calc_bond_forces_kernel(int start, int end, bond_t* bonds, coord_t* coords, cbond_t* cbonds, dvel_t* dvelocities, energy_accum_t* energy_sum) {
    int idx = blockIdx.x * blockDim.x + threadIdx.x + start;
    if (idx >= end) return;
    bond_t bond = bonds[idx];
    coord_t ri = coords[bond.ai - 1];
    coord_t rj = coords[bond.aj - 1];
    cbond_t cbond = cbonds[bond.code - 1];

    const double dx = rj.x - ri.x;
    const double dy = rj.y - ri.y;
    const double dz = rj.z - ri.z;
    const double r = sqrt(dx * dx + dy * dy + dz * dz);

    const double dr = r - cbond.b0;
    const double energy = 0.5 * cbond.kb * dr * dr;

    atomic_add_energy(energy_sum, energy);

    // update forces
    const double f = cbond.kb * dr / r;
    atomic_add_force(&dvelocities[bond.aj - 1].x, f * dx);
    atomic_add_force(&dvelocities[bond.aj - 1].y, f * dy);
    atomic_add_force(&dvelocities[bond.aj - 1].z, f * dz);
    atomic_add_force(&dvelocities[bond.ai - 1].x, -f * dx);
    atomic_add_force(&dvelocities[bond.ai - 1].y, -f * dy);
    atomic_add_force(&dvelocities[bond.ai - 1].z, -f * dz);
}

double calc_bond_forces_host(int start, int end) {
    int N = end - start;
    if (N <= 0) return 0.0;
    using namespace CudaBondForce;
    int blockSize = 256;
    int numBlocks = (N + blockSize - 1) / blockSize;

    energy_accum_t energy = 0;
    cudaMemcpy(d_energy_sum, &energy, sizeof(energy_accum_t), cudaMemcpyHostToDevice);

    auto& host_ctx = Context::instance();
    bond_t* d_bonds = host_ctx.bonds->gpu_data_p;
    coord_t* d_coords = host_ctx.coords->gpu_data_p;
    cbond_t* d_cbonds = host_ctx.cbonds->gpu_data_p;
    dvel_t* d_dvelocities = host_ctx.dvelocities->gpu_data_p;

    calc_bond_forces_kernel<<<numBlocks, blockSize>>>(start, end, d_bonds, d_coords, d_cbonds, d_dvelocities, d_energy_sum);
    cudaDeviceSynchronize();
    cudaMemcpy(&energy, d_energy_sum, sizeof(energy_accum_t), cudaMemcpyDeviceToHost);

    return energy_from_accum(energy);
}

void init_bond_force_kernel_data() {
    using namespace CudaBondForce;
    if (!is_initialized) {
        check_cudaMalloc((void**)&d_energy_sum, sizeof(energy_accum_t));
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
