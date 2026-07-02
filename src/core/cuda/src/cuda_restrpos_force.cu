#include <iostream>
#include <vector>

#include "common/include/context.h"
#include "cuda/include/cuda_restrpos_force.cuh"
#include "cuda/include/cuda_utility.cuh"
#include "cuda_force_accumulation.cuh"

namespace CudaRestrposForce {
}  // namespace CudaRestrposForce

__global__ void calc_restrpos_forces_kernel(
    restrpos_t* restrspos,
    int n_restrspos,
    coord_t* coords,
    double* lambdas,
    int n_lambdas,
    energy_accum_t* energy,
    dvel_t* dvelocities) {
    int idx = blockIdx.x * blockDim.x + threadIdx.x;
    if (idx >= n_restrspos) return;
    int ir = idx;

    int state, i;

    state = restrspos[ir].ipsi - 1;
    i = restrspos[ir].a - 1;

    const double dx = coords[i].x - restrspos[ir].x.x;
    const double dy = coords[i].y - restrspos[ir].x.y;
    const double dz = coords[i].z - restrspos[ir].x.z;

    double lambda;
    if (restrspos[ir].ipsi != 0) {
        lambda = lambdas[state];
    } else {
        lambda = 1.0;
    }

    const double kx = restrspos[ir].k.x;
    const double ky = restrspos[ir].k.y;
    const double kz = restrspos[ir].k.z;

    const double ener = 0.5 * kx * dx * dx + 0.5 * ky * dy * dy + 0.5 * kz * dz * dz;

    atomic_add_force(&dvelocities[i].x, kx * dx * lambda);
    atomic_add_force(&dvelocities[i].y, ky * dy * lambda);
    atomic_add_force(&dvelocities[i].z, kz * dz * lambda);

    if (restrspos[ir].ipsi == 0) {
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
void calc_restrpos_forces_host() {
    auto& host = Context::instance();
    if (host.n_restrspos == 0) return;
    using namespace CudaRestrposForce;

    auto d_restrspos = host.restrspos->gpu_data_p;
    auto d_coords = host.coords->gpu_data_p;
    auto d_lambdas = host.lambdas->gpu_data_p;
    auto d_dvelocities = host.dvelocities->gpu_data_p;

    int blockSize = 256;
    int numBlocks = (host.n_restrspos + blockSize - 1) / blockSize;
    calc_restrpos_forces_kernel<<<numBlocks, blockSize>>>(
        d_restrspos,
        host.n_restrspos,
        d_coords,
        d_lambdas,
        host.n_lambdas,
        host.energy.device(),
        d_dvelocities);
}

void init_restrpos_force_kernel_data() {
}

void cleanup_restrpos_force() {
}
