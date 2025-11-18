#include "cuda/include/CudaAngleForce.cuh"
#include "utils.h"
#include "cuda/include/CudaUtility.cuh"

namespace CudaAngleForce {
bool is_initialized = false;
angle_t* d_angles = nullptr;
coord_t* d_coords = nullptr;
cangle_t* d_cangles = nullptr;
dvel_t* d_dvelocities = nullptr;
double* d_energy_sum = nullptr;
}  // namespace CudaAngleForce


__global__ void calc_angle_forces_kernel(int start, int end, angle_t* angles, coord_t* coords, cangle_t* cangles, dvel_t* dvelocities, double* energy_sum) {
    int idx = blockIdx.x * blockDim.x + threadIdx.x + start;
    if (idx >= end) return;

    int i = angles[idx].ai - 1;
    int j = angles[idx].aj - 1;
    int k = angles[idx].ak - 1;

    coord_t ri = coords[i];
    coord_t rj = coords[j];
    coord_t rk = coords[k];

    cangle_t cang = cangles[angles[idx].code - 1];

    coord_t rji = {ri.x - rj.x, ri.y - rj.y, ri.z - rj.z};
    coord_t rjk = {rk.x - rj.x, rk.y - rj.y, rk.z - rj.z};

    double rji_length = sqrt(rji.x * rji.x + rji.y * rji.y + rji.z * rji.z);
    double rjk_length = sqrt(rjk.x * rjk.x + rjk.y * rjk.y + rjk.z * rjk.z);

    double cos_theta = (rji.x * rjk.x + rji.y * rjk.y + rji.z * rjk.z) / (rji_length * rjk_length);

    cos_theta = fmax(fmin(cos_theta, 1.0), -1.0);  // Clamp value to avoid NaNs
    double theta = acos(cos_theta);

    double dtheta = theta - to_radians_device(cang.th0);
    double energy = 0.5 * cang.kth * dtheta * dtheta;

    // calculate force magnitude
    double dv = cang.kth * dtheta;

    double f1 = sin(theta);
    if (fabs(f1) < 1e-12) {
        f1 = -1.0e12;
    } else {
        f1 = -1.0 / f1;
    }

    atomicAdd(energy_sum, energy);

    coord_t di = {
        f1 * (rjk.x / (rji_length * rjk_length) - cos_theta * rji.x / (rji_length * rji_length)),
        f1 * (rjk.y / (rji_length * rjk_length) - cos_theta * rji.y / (rji_length * rji_length)),
        f1 * (rjk.z / (rji_length * rjk_length) - cos_theta * rji.z / (rji_length * rji_length))};

    coord_t dk = {
        f1 * (rji.x / (rji_length * rjk_length) - cos_theta * rjk.x / (rjk_length * rjk_length)),
        f1 * (rji.y / (rji_length * rjk_length) - cos_theta * rjk.y / (rjk_length * rjk_length)),
        f1 * (rji.z / (rji_length * rjk_length) - cos_theta * rjk.z / (rjk_length * rjk_length))};

    atomicAdd(&dvelocities[i].x, dv * di.x);
    atomicAdd(&dvelocities[i].y, dv * di.y);
    atomicAdd(&dvelocities[i].z, dv * di.z);

    atomicAdd(&dvelocities[k].x, dv * dk.x);
    atomicAdd(&dvelocities[k].y, dv * dk.y);
    atomicAdd(&dvelocities[k].z, dv * dk.z);

    atomicAdd(&dvelocities[j].x, -dv * (di.x + dk.x));
    atomicAdd(&dvelocities[j].y, -dv * (di.y + dk.y));
    atomicAdd(&dvelocities[j].z, -dv * (di.z + dk.z));
}

double calc_angle_forces_host(int start, int end) {
    using namespace CudaAngleForce;
    int N = end - start;
    if (N <= 0) return 0.0;
    int blockSize = 256;
    int numBlocks = (N + blockSize - 1) / blockSize;

    if (!is_initialized) {
        is_initialized = true;
        check_cudaMalloc((void**)&d_angles, n_angles * sizeof(angle_t));
        check_cudaMalloc((void**)&d_coords, n_atoms * sizeof(coord_t));
        check_cudaMalloc((void**)&d_cangles, n_cangles * sizeof(cangle_t));
        check_cudaMalloc((void**)&d_dvelocities, n_atoms * sizeof(dvel_t));
        check_cudaMalloc((void**)&d_energy_sum, sizeof(double));
    }

    // copy data to device
    cudaMemcpy(d_angles, angles, n_angles * sizeof(angle_t), cudaMemcpyHostToDevice);
    cudaMemcpy(d_coords, coords, n_atoms * sizeof(coord_t), cudaMemcpyHostToDevice);
    cudaMemcpy(d_cangles, cangles, n_cangles * sizeof(cangle_t), cudaMemcpyHostToDevice);
    cudaMemcpy(d_dvelocities, dvelocities, n_atoms * sizeof(dvel_t), cudaMemcpyHostToDevice);

    double h_energy_sum = 0.0;
    cudaMemcpy(d_energy_sum, &h_energy_sum, sizeof(double), cudaMemcpyHostToDevice);

    // launch kernel
    calc_angle_forces_kernel<<<numBlocks, blockSize>>>(start, end, d_angles, d_coords, d_cangles, d_dvelocities, d_energy_sum);
    cudaDeviceSynchronize();

    // copy results back to host
    cudaMemcpy(&h_energy_sum, d_energy_sum, sizeof(double), cudaMemcpyDeviceToHost);
    cudaMemcpy(dvelocities, d_dvelocities, n_atoms * sizeof(dvel_t), cudaMemcpyDeviceToHost);
    return h_energy_sum;
}