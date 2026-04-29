#include "cuda/include/cuda_restrang_force.cuh"
#include "cuda/include/cuda_utility.cuh"
#include "common/include/context.h"
namespace CudaRestrangForce {
bool is_initialized = false;
real_t* d_E_restraint;
}  // namespace CudaRestrangForce

__global__ void calc_restrang_force_kernel(
    restrang_t* restrangs,
    int n_restrangs,
    coord_t* coords,
    real_t* lambdas,
    int n_lambdas,
    dvel_t* dvelocities,
    E_restraint_t* EQ_restraint,
    real_t* E_restraint) {
    int idx = blockIdx.x * blockDim.x + threadIdx.x;
    if (idx >= n_restrangs) return;
    int ir = idx;

    int state, i, j, k;
    coord_t dr, dr2, di, dk;
    real_t lambda, r2ij, r2jk, rij, rjk, cos_th, th;
    real_t dth, dv, ener, f1;

    state = restrangs[ir].ipsi - 1;
    i = restrangs[ir].ai - 1;
    j = restrangs[ir].aj - 1;
    k = restrangs[ir].ak - 1;

    // distance from atom i to atom j
    dr.x = coords[i].x - coords[j].x;
    dr.y = coords[i].y - coords[j].y;
    dr.z = coords[i].z - coords[j].z;

    // distance from atom k to atom j
    dr2.x = coords[k].x - coords[j].x;
    dr2.y = coords[k].y - coords[j].y;
    dr2.z = coords[k].z - coords[j].z;

    if (restrangs[ir].ipsi != 0) {
        lambda = lambdas[state];
    } else {
        lambda = 1;
    }

    r2ij = dr.x * dr.x + dr.y * dr.y + dr.z * dr.z;
    r2jk = dr2.x * dr2.x + dr2.y * dr2.y + dr2.z * dr2.z;

    rij = sqrt(r2ij);
    rjk = sqrt(r2jk);

    cos_th = dr.x * dr2.x + dr.y * dr2.y + dr.z * dr2.z;
    cos_th /= rij * rjk;

    if (cos_th > 1) cos_th = 1;
    if (cos_th < -1) cos_th = -1;

    th = acos(cos_th);
    dth = th - to_radians_device(restrangs[ir].ang);

    ener = .5 * restrangs[ir].k * dth * dth;
    dv = lambda * restrangs[ir].k * dth;

    f1 = sin(th);
    if (fabs(f1) < k_singular_sin_epsilon) {
        f1 = -1.0 / k_singular_sin_epsilon;
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

    if (restrangs[ir].ipsi == 0) {
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
    auto& host = Context::instance();
    if (host.n_restrangs == 0) return;
    using namespace CudaRestrangForce;

    auto d_restrangs = host.restrangs->gpu_data_p;
    auto d_coords = host.coords->gpu_data_p;
    auto d_lambdas = host.lambdas->gpu_data_p;
    auto d_dvelocities = host.dvelocities->gpu_data_p;
    auto d_EQ_restraint = host.EQ_restraint->gpu_data_p;

    real_t val = 0;
    cudaMemcpy(d_E_restraint, &val, sizeof(real_t), cudaMemcpyHostToDevice);

    int blockSize = 256;
    int numBlocks = (host.n_restrangs + blockSize - 1) / blockSize;
    calc_restrang_force_kernel<<<numBlocks, blockSize>>>(
        d_restrangs,
        host.n_restrangs,
        d_coords,
        d_lambdas,
        host.n_lambdas,
        d_dvelocities,
        d_EQ_restraint,
        d_E_restraint);
    cudaDeviceSynchronize();
    host.EQ_restraint->download();
    cudaMemcpy(&val, d_E_restraint, sizeof(real_t), cudaMemcpyDeviceToHost);
    host.E_restraint.Upres += val;
}

void init_restrang_force_kernel_data() {
    using namespace CudaRestrangForce;
    if (!is_initialized) {
        check_cudaMalloc((void**)&d_E_restraint, sizeof(real_t));
        is_initialized = true;
    }
}

void cleanup_restrang_force() {
    using namespace CudaRestrangForce;
    if (is_initialized) {
        cudaFree(d_E_restraint);
        is_initialized = false;
    }
}
