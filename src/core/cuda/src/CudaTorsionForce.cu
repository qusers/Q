#include "cuda/include/CudaTorsionForce.cuh"
#include "cuda/include/CudaUtility.cuh"
#include "utils.h"

namespace CudaTorsionForce {
bool is_initialized = false;
torsion_t* d_torsions = nullptr;
ctorsion_t* d_ctorsions = nullptr;
coord_t* d_coords = nullptr;
dvel_t* d_dvelocities = nullptr;
double* d_energy_sum = nullptr;
}  // namespace CudaTorsionForce

__global__ void calc_torsion_forces_kernel(int start, int end, torsion_t* torsions, ctorsion_t* ctorsions, coord_t* coords, dvel_t* dvelocities, double* energy_sum) {
    int i = blockIdx.x * blockDim.x + threadIdx.x + start;
    if (i >= end) return;
    int aii, aji, aki, ali;

    coord_t ai, aj, ak, al;
    coord_t rji, rjk, rkl, rnj, rnk, rki, rlj;
    coord_t di, dl, dpi, dpj, dpk, dpl;

    double bj2inv, bk2inv, bjinv, bkinv;
    double cos_phi, phi;
    double arg, dv, f1;
    double ener;
    double torsion = 0;

    torsion_t t;
    ctorsion_t ctors;

    t = torsions[i];
    ctors = ctorsions[t.code - 1];

    aii = t.ai - 1;
    aji = t.aj - 1;
    aki = t.ak - 1;
    ali = t.al - 1;

    ai = coords[aii];
    aj = coords[aji];
    ak = coords[aki];
    al = coords[ali];

    rji.x = ai.x - aj.x;
    rji.y = ai.y - aj.y;
    rji.z = ai.z - aj.z;

    rjk.x = ak.x - aj.x;
    rjk.y = ak.y - aj.y;
    rjk.z = ak.z - aj.z;

    rkl.x = al.x - ak.x;
    rkl.y = al.y - ak.y;
    rkl.z = al.z - ak.z;

    rnj.x = rji.y * rjk.z - rji.z * rjk.y;
    rnj.y = rji.z * rjk.x - rji.x * rjk.z;
    rnj.z = rji.x * rjk.y - rji.y * rjk.x;

    rnk.x = -rjk.y * rkl.z + rjk.z * rkl.y;
    rnk.y = -rjk.z * rkl.x + rjk.x * rkl.z;
    rnk.z = -rjk.x * rkl.y + rjk.y * rkl.x;

    bj2inv = 1 / (pow(rnj.x, 2) + pow(rnj.y, 2) + pow(rnj.z, 2));
    bk2inv = 1 / (pow(rnk.x, 2) + pow(rnk.y, 2) + pow(rnk.z, 2));
    bjinv = sqrt(bj2inv);
    bkinv = sqrt(bk2inv);

    cos_phi = (rnj.x * rnk.x + rnj.y * rnk.y + rnj.z * rnk.z) * (bjinv * bkinv);
    cos_phi = fmin(fmax(cos_phi, -1.0), 1.0);
    phi = acos(cos_phi);
    if (rjk.x * (rnj.y * rnk.z - rnj.z * rnk.y) + rjk.y * (rnj.z * rnk.x - rnj.x * rnk.z) + rjk.z * (rnj.x * rnk.y - rnj.y * rnk.x) < 0) {
        phi = -phi;
    }

    // Energy
    arg = ctors.n * phi - to_radians_device(ctors.d);
    ener = ctors.k * (1 + cos(arg));
    dv = -ctors.n * ctors.k * sin(arg);

    // Forces
    f1 = sin(phi);
    if (abs(f1) < 1E-12) f1 = 1E-12;
    f1 = -1 / f1;

    di.x = f1 * (rnk.x * (bjinv * bkinv) - cos_phi * rnj.x * bj2inv);
    di.y = f1 * (rnk.y * (bjinv * bkinv) - cos_phi * rnj.y * bj2inv);
    di.z = f1 * (rnk.z * (bjinv * bkinv) - cos_phi * rnj.z * bj2inv);
    dl.x = f1 * (rnj.x * (bjinv * bkinv) - cos_phi * rnk.x * bk2inv);
    dl.y = f1 * (rnj.y * (bjinv * bkinv) - cos_phi * rnk.y * bk2inv);
    dl.z = f1 * (rnj.z * (bjinv * bkinv) - cos_phi * rnk.z * bk2inv);

    rki.x = rji.x - rjk.x;
    rki.y = rji.y - rjk.y;
    rki.z = rji.z - rjk.z;
    rlj.x = -rjk.x - rkl.x;
    rlj.y = -rjk.y - rkl.y;
    rlj.z = -rjk.z - rkl.z;

    dpi.x = rjk.y * di.z - rjk.z * di.y;
    dpi.y = rjk.z * di.x - rjk.x * di.z;
    dpi.z = rjk.x * di.y - rjk.y * di.x;
    dpj.x = rki.y * di.z - rki.z * di.y + rkl.y * dl.z - rkl.z * dl.y;
    dpj.y = rki.z * di.x - rki.x * di.z + rkl.z * dl.x - rkl.x * dl.z;
    dpj.z = rki.x * di.y - rki.y * di.x + rkl.x * dl.y - rkl.y * dl.x;
    dpk.x = rlj.y * dl.z - rlj.z * dl.y - rji.y * di.z + rji.z * di.y;
    dpk.y = rlj.z * dl.x - rlj.x * dl.z - rji.z * di.x + rji.x * di.z;
    dpk.z = rlj.x * dl.y - rlj.y * dl.x - rji.x * di.y + rji.y * di.x;
    dpl.x = rjk.y * dl.z - rjk.z * dl.y;
    dpl.y = rjk.z * dl.x - rjk.x * dl.z;
    dpl.z = rjk.x * dl.y - rjk.y * dl.x;

    // Update energy and forces
    atomicAdd(energy_sum, ener);

    atomicAdd(&dvelocities[aii].x, dv * dpi.x);
    atomicAdd(&dvelocities[aii].y, dv * dpi.y);
    atomicAdd(&dvelocities[aii].z, dv * dpi.z);
    atomicAdd(&dvelocities[aji].x, dv * dpj.x);
    atomicAdd(&dvelocities[aji].y, dv * dpj.y);
    atomicAdd(&dvelocities[aji].z, dv * dpj.z);
    atomicAdd(&dvelocities[aki].x, dv * dpk.x);
    atomicAdd(&dvelocities[aki].y, dv * dpk.y);
    atomicAdd(&dvelocities[aki].z, dv * dpk.z);
    atomicAdd(&dvelocities[ali].x, dv * dpl.x);
    atomicAdd(&dvelocities[ali].y, dv * dpl.y);
    atomicAdd(&dvelocities[ali].z, dv * dpl.z);
}

double calc_torsion_forces_host(int start, int end) {
    using namespace CudaTorsionForce;
    int N = end - start;
    if (N <= 0) return 0.0;
    int blockSize = 256;
    int numBlocks = (N + blockSize - 1) / blockSize;
    if (!is_initialized) {
        check_cudaMalloc((void**)&d_torsions, sizeof(torsion_t) * n_torsions);
        check_cudaMalloc((void**)&d_ctorsions, sizeof(ctorsion_t) * n_ctorsions);
        check_cudaMalloc((void**)&d_coords, sizeof(coord_t) * n_atoms);
        check_cudaMalloc((void**)&d_dvelocities, sizeof(dvel_t) * n_atoms);
        check_cudaMalloc((void**)&d_energy_sum, sizeof(double));
        is_initialized = true;
    }

    cudaMemcpy(d_torsions, torsions, sizeof(torsion_t) * n_torsions, cudaMemcpyHostToDevice);
    cudaMemcpy(d_ctorsions, ctorsions, sizeof(ctorsion_t) * n_ctorsions, cudaMemcpyHostToDevice);
    cudaMemcpy(d_coords, coords, sizeof(coord_t) * n_atoms, cudaMemcpyHostToDevice);
    cudaMemcpy(d_dvelocities, dvelocities, sizeof(dvel_t) * n_atoms, cudaMemcpyHostToDevice);
    double zero = 0.0;
    cudaMemcpy(d_energy_sum, &zero, sizeof(double), cudaMemcpyHostToDevice);
    calc_torsion_forces_kernel<<<numBlocks, blockSize>>>(start, end, d_torsions, d_ctorsions, d_coords, d_dvelocities, d_energy_sum);
    cudaDeviceSynchronize();
    cudaMemcpy(&zero, d_energy_sum, sizeof(double), cudaMemcpyDeviceToHost);
    cudaMemcpy(dvelocities, d_dvelocities, sizeof(dvel_t) * n_atoms, cudaMemcpyDeviceToHost);
    return zero;
}
