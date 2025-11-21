#include "cuda/include/CudaRestrposForce.cuh"
#include "cuda/include/CudaUtility.cuh"
#include "utils.h"

namespace CudaRestrposForce {
bool is_initialized = false;
restrpos_t* d_restrpos = nullptr;
coord_t* d_coords = nullptr;
double* d_lambda = nullptr;
E_restraint_t* d_EQ_restraint = nullptr;
double* d_E_restraint = nullptr;
dvel_t* d_dvelocities = nullptr;

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
    using namespace CudaRestrposForce;
    if (!is_initialized) {
        check_cudaMalloc((void**)&d_restrpos, sizeof(restrpos_t) * n_restrspos);
        check_cudaMalloc((void**)&d_coords, sizeof(coord_t) * n_atoms);
        check_cudaMalloc((void**)&d_lambda, sizeof(double));
        check_cudaMalloc((void**)&d_EQ_restraint, sizeof(E_restraint_t) * n_lambdas);
        check_cudaMalloc((void**)&d_E_restraint, sizeof(double));
        check_cudaMalloc((void**)&d_dvelocities, sizeof(dvel_t) * n_atoms);

        cudaMemcpy(d_restrpos, restrspos, sizeof(restrpos_t) * n_restrspos, cudaMemcpyHostToDevice);
        cudaMemcpy(d_lambda, &lambdas, sizeof(double), cudaMemcpyHostToDevice);
        is_initialized = true;
    }
    cudaMemcpy(d_coords, coords, sizeof(coord_t) * n_atoms, cudaMemcpyHostToDevice);
    cudaMemcpy(d_dvelocities, dvelocities, sizeof(dvel_t) * n_atoms, cudaMemcpyHostToDevice);
    cudaMemcpy(d_EQ_restraint, EQ_restraint, sizeof(E_restraint_t) * n_lambdas, cudaMemcpyHostToDevice);
    double val = 0.0;
    cudaMemcpy(d_E_restraint, &val, sizeof(double), cudaMemcpyHostToDevice);

    int blockSize = 256;
    int numBlocks = (n_restrspos + blockSize - 1) / blockSize;
    calc_restrpos_forces_kernel<<<numBlocks, blockSize>>>(
        d_restrpos,
        n_restrspos,
        d_coords,
        d_lambda,
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
void cleanup_restrpos_force() {
    using namespace CudaRestrposForce;
    if (is_initialized) {
        cudaFree(d_restrpos);
        cudaFree(d_coords);
        cudaFree(d_lambda);
        cudaFree(d_EQ_restraint);
        cudaFree(d_E_restraint);
        cudaFree(d_dvelocities);
        is_initialized = false;
    }
}
