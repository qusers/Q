#include <iostream>

#include "cuda/include/cuda_restrpos_force.cuh"
#include "cuda/include/cuda_utility.cuh"
#include "common/include/context.h"
#include "cuda_force_accumulation.cuh"

namespace CudaRestrposForce {
bool is_initialized = false;
real_t* d_E_restraint;
}  // namespace CudaRestrposForce

__global__ void calc_restrpos_forces_kernel(
    restrpos_t* restrspos,
    int n_restrspos,
    coord_t* coords,
    real_t* lambdas,
    int n_lambdas,
    E_restraint_t* EQ_restraint,
    real_t* E_restraint,
    dvel_t* dvelocities) {
    int idx = blockIdx.x * blockDim.x + threadIdx.x;
    if (idx >= n_restrspos) return;
    int ir = idx;

    int state, i;
    coord_t dr;
    real_t lambda, ener, x2, y2, z2;

    state = restrspos[ir].ipsi - 1;
    i = restrspos[ir].a - 1;

    dr.x = coords[i].x - restrspos[ir].x.x;
    dr.y = coords[i].y - restrspos[ir].x.y;
    dr.z = coords[i].z - restrspos[ir].x.z;

    if (restrspos[ir].ipsi != 0) {
        lambda = lambdas[state];
    } else {
        lambda = 1;
    }

    x2 = dr.x * dr.x;
    y2 = dr.y * dr.y;
    z2 = dr.z * dr.z;

    ener = .5 * restrspos[ir].k.x * x2 + .5 * restrspos[ir].k.y * y2 + .5 * restrspos[ir].k.z * z2;

    atomic_add_force(&dvelocities[i].x, restrspos[ir].k.x * dr.x * lambda);
    atomic_add_force(&dvelocities[i].y, restrspos[ir].k.y * dr.y * lambda);
    atomic_add_force(&dvelocities[i].z, restrspos[ir].k.z * dr.z * lambda);

    if (restrspos[ir].ipsi == 0) {
        for (int k = 0; k < n_lambdas; k++) {
            atomicAdd(&EQ_restraint[k].Urestr, ener);
        }
        if (n_lambdas == 0) {
            atomicAdd(E_restraint, ener);
        }
    } else {
        atomicAdd(&EQ_restraint[state].Urestr, ener);
    }
}
void calc_restrpos_forces_host() {
    auto& host = Context::instance();
    if (host.n_restrspos == 0) return;
    using namespace CudaRestrposForce;
    real_t val = 0.0;
    cudaMemcpy(d_E_restraint, &val, sizeof(real_t), cudaMemcpyHostToDevice);

    auto d_restrspos = host.restrspos->gpu_data_p;
    auto d_coords = host.coords->gpu_data_p;
    auto d_lambdas = host.lambdas->gpu_data_p;
    auto d_EQ_restraint = host.EQ_restraint->gpu_data_p;
    auto d_dvelocities = host.dvelocities->gpu_data_p;

    int blockSize = 256;
    int numBlocks = (host.n_restrspos + blockSize - 1) / blockSize;
    calc_restrpos_forces_kernel<<<numBlocks, blockSize>>>(
        d_restrspos,
        host.n_restrspos,
        d_coords,
        d_lambdas,
        host.n_lambdas,
        d_EQ_restraint,
        d_E_restraint,
        d_dvelocities);
    cudaDeviceSynchronize();
    cudaMemcpy(&val, d_E_restraint, sizeof(real_t), cudaMemcpyDeviceToHost);
    host.E_restraint.Upres += val;
    host.EQ_restraint->download();
}

void init_restrpos_force_kernel_data() {
    using namespace CudaRestrposForce;
    if (!is_initialized) {
        check_cudaMalloc((void**)&d_E_restraint, sizeof(real_t));
        is_initialized = true;
    }
}

void cleanup_restrpos_force() {
    using namespace CudaRestrposForce;
    if (is_initialized) {
        cudaFree(d_E_restraint);
        is_initialized = false;
    }
}
