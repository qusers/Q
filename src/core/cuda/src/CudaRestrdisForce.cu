#include "cuda/include/CudaRestrdisForce.cuh"
#include "cuda/include/CudaUtility.cuh"
#include "utils.h"
#include <iostream>

namespace CudaRestrdisForce {
bool is_initialized = false;
restrdis_t* d_restrdists = nullptr;
coord_t* d_coords = nullptr;
double* d_lambdas = nullptr;
dvel_t* d_dvelocities = nullptr;
E_restraint_t* d_EQ_restraint = nullptr;
double* d_E_restraint = nullptr;

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
    using namespace CudaRestrdisForce;
    if (!is_initialized) {
        check_cudaMalloc((void**)&d_restrdists, n_restrdists * sizeof(restrdis_t));
        check_cudaMalloc((void**)&d_coords, n_atoms * sizeof(coord_t));
        check_cudaMalloc((void**)&d_lambdas, n_lambdas * sizeof(double));
        check_cudaMalloc((void**)&d_dvelocities, n_atoms * sizeof(dvel_t));
        check_cudaMalloc((void**)&d_EQ_restraint, sizeof(E_restraint_t) * n_lambdas);
        check_cudaMalloc((void**)&d_E_restraint, sizeof(double) * n_lambdas);

        cudaMemcpy(d_restrdists, restrdists, n_restrdists * sizeof(restrdis_t), cudaMemcpyHostToDevice);
        cudaMemcpy(d_lambdas, lambdas, n_lambdas * sizeof(double), cudaMemcpyHostToDevice);
        is_initialized = true;
    }
    cudaMemcpy(d_coords, coords, n_atoms * sizeof(coord_t), cudaMemcpyHostToDevice);
    cudaMemcpy(d_dvelocities, dvelocities, n_atoms * sizeof(dvel_t), cudaMemcpyHostToDevice);
    cudaMemcpy(d_EQ_restraint, EQ_restraint, sizeof(E_restraint_t) * n_lambdas, cudaMemcpyHostToDevice);
    cudaMemset(d_E_restraint, 0, sizeof(double) * n_lambdas);
    int blockSize = 256;
    int numBlocks = (n_restrdists + blockSize - 1) / blockSize;
    calc_restrdis_forces_kernel<<<numBlocks, blockSize>>>(
        d_restrdists,
        n_restrdists,
        d_coords,
        d_lambdas,
        n_lambdas,
        d_dvelocities,
        d_EQ_restraint,
        d_E_restraint);
    cudaDeviceSynchronize();
    cudaMemcpy(dvelocities, d_dvelocities, n_atoms * sizeof(dvel_t), cudaMemcpyDeviceToHost);
    cudaMemcpy(EQ_restraint, d_EQ_restraint, sizeof(E_restraint_t) * n_lambdas, cudaMemcpyDeviceToHost);
    double ener;
    cudaMemcpy(&ener, d_E_restraint, sizeof(double), cudaMemcpyDeviceToHost);
    printf("Energy restraint: %f\n", ener);
    E_restraint.Urestr += ener;
}

void cleanup_restrdis_force() {
    using namespace CudaRestrdisForce;
    if (is_initialized) {
        cudaFree(d_restrdists);
        cudaFree(d_coords);
        cudaFree(d_lambdas);
        cudaFree(d_dvelocities);
        cudaFree(d_EQ_restraint);
        cudaFree(d_E_restraint);
        is_initialized = false;
    }
}
