#include "cuda/include/CudaContext.cuh"
#include "cuda/include/CudaRestrangForce.cuh"
#include "cuda/include/CudaUtility.cuh"
#include "utils.h"

__global__ void calc_restrang_force_kernel(
    restrang_t* restrangs,
    int n_restrangs,
    coord_t* coords,
    int n_atoms,
    double* lambdas,
    int n_lambdas,
    dvel_t* dvelocities,
    E_restraint_t* EQ_restraint,
    double* E_restraint,
    restrdis_t* restrdists) {
    int idx = blockIdx.x * blockDim.x + threadIdx.x;
    if (idx >= n_restrangs) return;
    int ir = idx;

    int state, i, j, k;
    coord_t dr, dr2, di, dk;
    double lambda, r2ij, r2jk, rij, rjk, cos_th, th;
    double dth, dv, ener, f1;

    state = restrangs[ir].ipsi - 1;
    i = restrangs[ir].ai - 1;
    j = restrangs[ir].aj - 1;
    k = restrangs[ir].ak - 1;

    // distance from atom i to atom j
    dr.x = coords[i].x - coords[j].x;
    dr.y = coords[i].y - coords[j].y;
    dr.z = coords[i].z - coords[j].z;

    // distance from atom k to atom j
    dr.x = coords[k].x - coords[j].x;
    dr.y = coords[k].y - coords[j].y;
    dr.z = coords[k].z - coords[j].z;

    if (restrangs[ir].ipsi != 0) {
        lambda = lambdas[state];
    } else {
        lambda = 1;
    }

    r2ij = pow(dr.x, 2) + pow(dr.y, 2) + pow(dr.z, 2);
    r2jk = pow(dr2.x, 2) + pow(dr2.y, 2) + pow(dr2.z, 2);

    rij = sqrt(r2ij);
    rjk = sqrt(r2jk);

    cos_th = dr.x * dr2.x + dr.y * dr2.y + dr.z * dr2.z;
    cos_th /= rij * rjk;

    if (cos_th > 1) cos_th = 1;
    if (cos_th < -1) cos_th = -1;

    th = acos(cos_th);
    dth = th - to_radians_device(restrangs[ir].ang);

    ener = .5 * restrangs[ir].k * pow(dth, 2);
    dv = lambda * restrangs[ir].k * dth;

    f1 = sin(th);
    if (abs(f1) < 1E-12) {
        f1 = -1E-12;
    } else {
        f1 = -1 / f1;
    }

    di.x = f1 * (dr2.x / (rij * rjk) - cos_th * dr.x / r2ij);
    di.y = f1 * (dr2.y / (rij * rjk) - cos_th * dr.y / r2ij);
    di.z = f1 * (dr2.z / (rij * rjk) - cos_th * dr.z / r2ij);
    dk.x = f1 * (dr.x / (rij * rjk) - cos_th * dr2.x / r2jk);
    dk.y = f1 * (dr.y / (rij * rjk) - cos_th * dr2.y / r2jk);
    dk.z = f1 * (dr.z / (rij * rjk) - cos_th * dr2.z / r2jk);

    atomicAdd(&dvelocities[i].x, dv * di.x);
    atomicAdd(&dvelocities[i].y, dv * di.y);
    atomicAdd(&dvelocities[i].z, dv * di.z);
    atomicAdd(&dvelocities[k].x, dv * dk.x);
    atomicAdd(&dvelocities[k].y, dv * dk.y);
    atomicAdd(&dvelocities[k].z, dv * dk.z);
    atomicAdd(&dvelocities[j].x, -dv * (di.x + dk.x));
    atomicAdd(&dvelocities[j].y, -dv * (di.y + dk.y));
    atomicAdd(&dvelocities[j].z, -dv * (di.z + dk.z));

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

void calc_restrang_force_host() {
    CudaContext& ctx = CudaContext::instance();

    auto d_restrangs = ctx.d_restrangs;
    auto d_coords = ctx.d_coords;
    auto d_lambdas = ctx.d_lambdas;
    auto d_dvelocities = ctx.d_dvelocities;
    auto d_EQ_restraint = ctx.d_EQ_restraint;
    auto d_restrdists = ctx.d_restrdists;
    ctx.sync_all_to_device();

    double* d_E_restraint;
    double val = 0;
    check_cudaMalloc((void**)&d_E_restraint, sizeof(double));
    cudaMemcpy(d_E_restraint, &val, sizeof(double), cudaMemcpyHostToDevice);

    int blockSize = 256;
    int numBlocks = (n_restrangs + blockSize - 1) / blockSize;
    calc_restrang_force_kernel<<<numBlocks, blockSize>>>(
        d_restrangs,
        n_restrangs,
        d_coords,
        n_atoms,
        d_lambdas,
        n_lambdas,
        d_dvelocities,
        d_EQ_restraint,
        d_E_restraint,
        d_restrdists);
    cudaDeviceSynchronize();
    cudaMemcpy(coords, d_coords, sizeof(coord_t) * n_atoms, cudaMemcpyDeviceToHost);
    cudaMemcpy(dvelocities, d_dvelocities, sizeof(dvel_t) * n_atoms, cudaMemcpyDeviceToHost);
    cudaMemcpy(EQ_restraint, d_EQ_restraint, sizeof(E_restraint_t) * n_lambdas, cudaMemcpyDeviceToHost);
    cudaMemcpy(&val, d_E_restraint, sizeof(double), cudaMemcpyDeviceToHost);
    E_restraint.Urestr += val;
}
