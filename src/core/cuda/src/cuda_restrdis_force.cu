#include <iostream>

#include "cuda/include/cuda_context.cuh"
#include "cuda/include/cuda_restrdis_force.cuh"
#include "cuda/include/cuda_utility.cuh"
#include "common/include/context.h"
namespace CudaRestrdisForce {
bool is_initialized = false;
double* d_E_restraint;
}  // namespace CudaRestrdisForce

__global__ void calc_restrdis_forces_kernel(
    restrdis_t* restrdists,
    int n_restrdists,
    coord_t* coords,
    double* lambdas,
    int n_lambdas,
    dvel_t* dvelocities,
    E_restraint_t* EQ_restraint,
    double* E_restraint) {
    int idx = blockIdx.x * blockDim.x + threadIdx.x;
    if (idx >= n_restrdists) return;

    int state, i, j;
    coord_t dr;
    double lambda, b, db, dv, ener;

    int ir = idx;

    state = restrdists[ir].ipsi - 1;
    i = restrdists[ir].ai - 1;
    j = restrdists[ir].aj - 1;

    dr.x = coords[j].x - coords[i].x;
    dr.y = coords[j].y - coords[i].y;
    dr.z = coords[j].z - coords[i].z;

    if (restrdists[ir].ipsi != 0) {
        lambda = lambdas[state];
    } else {
        lambda = 1;
    }

    b = sqrt(pow(dr.x, 2) + pow(dr.y, 2) + pow(dr.z, 2));
    if (b < restrdists[ir].d1) {
        db = b - restrdists[ir].d1;
    } else if (b > restrdists[ir].d2) {
        db = b - restrdists[ir].d2;
    } else {
        db = 0;
        return;
    }

    ener = .5 * restrdists[ir].k * pow(db, 2);
    dv = lambda * restrdists[ir].k * db / b;

    atomicAdd(&dvelocities[j].x, dr.x * dv);
    atomicAdd(&dvelocities[j].y, dr.y * dv);
    atomicAdd(&dvelocities[j].z, dr.z * dv);
    atomicAdd(&dvelocities[i].x, -dr.x * dv);
    atomicAdd(&dvelocities[i].y, -dr.y * dv);
    atomicAdd(&dvelocities[i].z, -dr.z * dv);

    if (restrdists[ir].ipsi == 0) {
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

void calc_restrdis_forces_host() {
    auto& host = Context::instance();
    if (host.n_restrdists == 0) return;
    using namespace CudaRestrdisForce;
    CudaContext& ctx = CudaContext::instance();
    auto d_restrdists = host.restrdists->gpu_data_p;
    auto d_coords = host.coords->gpu_data_p;
    auto d_lambdas = host.lambdas->gpu_data_p;
    auto d_dvelocities = host.dvelocities->gpu_data_p;
    auto d_EQ_restraint = host.EQ_restraint->gpu_data_p;

    cudaMemset(d_E_restraint, 0, sizeof(double));

    int blockSize = 256;
    int numBlocks = (host.n_restrdists + blockSize - 1) / blockSize;
    calc_restrdis_forces_kernel<<<numBlocks, blockSize>>>(
        d_restrdists,
        host.n_restrdists,
        d_coords,
        d_lambdas,
        host.n_lambdas,
        d_dvelocities,
        d_EQ_restraint,
        d_E_restraint);
    cudaDeviceSynchronize();
    host.EQ_restraint->download();
    double ener;
    cudaMemcpy(&ener, d_E_restraint, sizeof(double), cudaMemcpyDeviceToHost);
    printf("Energy restraint: %f\n", ener);
    host.E_restraint.Upres += ener;
}

void init_restrdis_force_kernel_data() {
    using namespace CudaRestrdisForce;
    if (!is_initialized) {
        check_cudaMalloc((void**)&d_E_restraint, sizeof(double));
        is_initialized = true;
    }
}

void cleanup_restrdis_force() {
    using namespace CudaRestrdisForce;
    if (is_initialized) {
        cudaFree(d_E_restraint);
        is_initialized = false;
    }
}
