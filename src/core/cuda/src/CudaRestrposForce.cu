#include <iostream>

#include "cuda/include/CudaContext.cuh"
#include "cuda/include/CudaRestrposForce.cuh"
#include "cuda/include/CudaUtility.cuh"
#include "utils.h"

namespace CudaRestrposForce {
bool is_initialized = false;
double* d_E_restraint;
}  // namespace CudaRestrposForce

__global__ void calc_restrpos_forces_kernel(
    restrpos_t* restrspos,
    int n_restrspos,
    coord_t* coords,
    double* lambdas,
    int n_lambdas,
    E_restraint_t* EQ_restraint,
    double* E_restraint,
    dvel_t* dvelocities) {
    int idx = blockIdx.x * blockDim.x + threadIdx.x;
    if (idx >= n_restrspos) return;
    int ir = idx;

    int state, i;
    coord_t dr;
    double lambda, ener, x2, y2, z2;

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

    x2 = pow(dr.x, 2);
    y2 = pow(dr.y, 2);
    z2 = pow(dr.z, 2);

    ener = .5 * restrspos[ir].k.x * x2 + .5 * restrspos[ir].k.y * y2 + .5 * restrspos[ir].k.z * z2;

    atomicAdd(&dvelocities[i].x, restrspos[ir].k.x * dr.x * lambda);
    atomicAdd(&dvelocities[i].y, restrspos[ir].k.y * dr.y * lambda);
    atomicAdd(&dvelocities[i].z, restrspos[ir].k.z * dr.z * lambda);

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
    if (n_restrspos == 0) return;
    using namespace CudaRestrposForce;
    double val = 0.0;
    cudaMemcpy(d_E_restraint, &val, sizeof(double), cudaMemcpyHostToDevice);

    CudaContext& ctx = CudaContext::instance();
    auto d_restrspos = ctx.d_restrspos;
    auto d_coords = ctx.d_coords;
    auto d_lambdas = ctx.d_lambdas;
    auto d_EQ_restraint = ctx.d_EQ_restraint;
    auto d_dvelocities = ctx.d_dvelocities;

    int blockSize = 256;
    int numBlocks = (n_restrspos + blockSize - 1) / blockSize;
    calc_restrpos_forces_kernel<<<numBlocks, blockSize>>>(
        d_restrspos,
        n_restrspos,
        d_coords,
        d_lambdas,
        n_lambdas,
        d_EQ_restraint,
        d_E_restraint,
        d_dvelocities);
    cudaDeviceSynchronize();
    cudaMemcpy(dvelocities, d_dvelocities, sizeof(dvel_t) * n_atoms, cudaMemcpyDeviceToHost);
    cudaMemcpy(&val, d_E_restraint, sizeof(double), cudaMemcpyDeviceToHost);
    E_restraint.Urestr += val;
    cudaMemcpy(EQ_restraint, d_EQ_restraint, sizeof(E_restraint_t) * n_lambdas, cudaMemcpyDeviceToHost);
}

void init_restrpos_force_kernel_data() {
    using namespace CudaRestrposForce;
    if (!is_initialized) {
        check_cudaMalloc((void**)&d_E_restraint, sizeof(double));
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