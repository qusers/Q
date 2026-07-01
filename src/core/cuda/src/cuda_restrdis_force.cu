#include <iostream>
#include <vector>

#include "common/include/context.h"
#include "cuda/include/cuda_restrdis_force.cuh"
#include "cuda/include/cuda_utility.cuh"
#include "cuda_force_accumulation.cuh"
namespace CudaRestrdisForce {
}  // namespace CudaRestrdisForce

__global__ void calc_restrdis_forces_kernel(
    restrdis_t* restrdists,
    int n_restrdists,
    coord_t* coords,
    double* lambdas,
    int n_lambdas,
    dvel_t* dvelocities,
    energy_accum_t* energy) {
    int idx = blockIdx.x * blockDim.x + threadIdx.x;
    if (idx >= n_restrdists) return;

    int state, i, j;

    int ir = idx;

    state = restrdists[ir].ipsi - 1;
    i = restrdists[ir].ai - 1;
    j = restrdists[ir].aj - 1;

    const double dx = coords[j].x - coords[i].x;
    const double dy = coords[j].y - coords[i].y;
    const double dz = coords[j].z - coords[i].z;

    double lambda;
    if (restrdists[ir].ipsi != 0) {
        lambda = lambdas[state];
    } else {
        lambda = 1.0;
    }

    const double b = sqrt(dx * dx + dy * dy + dz * dz);
    double db;
    if (b < restrdists[ir].d1) {
        db = b - restrdists[ir].d1;
    } else if (b > restrdists[ir].d2) {
        db = b - restrdists[ir].d2;
    } else {
        return;
    }

    const double ener = 0.5 * restrdists[ir].k * db * db;
    const double dv = lambda * restrdists[ir].k * db / b;

    atomic_add_force(&dvelocities[j].x, dx * dv);
    atomic_add_force(&dvelocities[j].y, dy * dv);
    atomic_add_force(&dvelocities[j].z, dz * dv);
    atomic_add_force(&dvelocities[i].x, -dx * dv);
    atomic_add_force(&dvelocities[i].y, -dy * dv);
    atomic_add_force(&dvelocities[i].z, -dz * dv);

    if (restrdists[ir].ipsi == 0) {
        for (int k = 0; k < n_lambdas; k++) {
            atomic_add_energy(&energy[EnergyBuffer::eq_index(ENERGY_FIXED_COUNT, k, EQ_RESTR_URESTR)], ener);
        }
        if (n_lambdas == 0) {
            atomic_add_energy(&energy[E_RESTR_PRES], ener);
        }
    } else {
        atomic_add_energy(&energy[EnergyBuffer::eq_index(ENERGY_FIXED_COUNT, state, EQ_RESTR_URESTR)], ener);
    }
}

void calc_restrdis_forces_host() {
    auto& host = Context::instance();
    if (host.n_restrdists == 0) return;
    using namespace CudaRestrdisForce;
    auto d_restrdists = host.restrdists->gpu_data_p;
    auto d_coords = host.coords->gpu_data_p;
    auto d_lambdas = host.lambdas->gpu_data_p;
    auto d_dvelocities = host.dvelocities->gpu_data_p;

    int blockSize = 256;
    int numBlocks = (host.n_restrdists + blockSize - 1) / blockSize;
    calc_restrdis_forces_kernel<<<numBlocks, blockSize>>>(
        d_restrdists,
        host.n_restrdists,
        d_coords,
        d_lambdas,
        host.n_lambdas,
        d_dvelocities,
        host.energy.device());
}

void init_restrdis_force_kernel_data() {
}

void cleanup_restrdis_force() {
}
