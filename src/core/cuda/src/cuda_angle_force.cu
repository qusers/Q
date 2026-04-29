#include "cuda/include/cuda_angle_force.cuh"
#include "cuda/include/cuda_utility.cuh"
#include "context.h"

namespace CudaAngleForce {
bool is_initialized = false;
real_t* d_energy_sum;
}  // namespace CudaAngleForce

__global__ void calc_angle_forces_kernel(int start, int end, angle_t* angles, coord_t* coords, cangle_t* cangles, dvel_t* dvelocities, real_t* energy_sum) {
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

    real_t rji_length = sqrt(rji.x * rji.x + rji.y * rji.y + rji.z * rji.z);
    real_t rjk_length = sqrt(rjk.x * rjk.x + rjk.y * rjk.y + rjk.z * rjk.z);

    real_t cos_theta = (rji.x * rjk.x + rji.y * rjk.y + rji.z * rjk.z) / (rji_length * rjk_length);

    cos_theta = cos_theta > static_cast<real_t>(1.0) ? static_cast<real_t>(1.0) : cos_theta;
    cos_theta = cos_theta < static_cast<real_t>(-1.0) ? static_cast<real_t>(-1.0) : cos_theta;
    real_t theta = acos(cos_theta);

    real_t dtheta = theta - to_radians_device(cang.th0);
    real_t energy = 0.5 * cang.kth * dtheta * dtheta;

    // calculate force magnitude
    real_t dv = cang.kth * dtheta;

    real_t f1 = sin(theta);
    if (fabs(f1) < k_singular_sin_epsilon) {
        f1 = -1.0 / k_singular_sin_epsilon;
    } else {
        f1 = -1.0 / f1;
    }

    atomicAdd(energy_sum, energy);

    coord_t di = {
        static_cast<real_t>(f1 * (rjk.x / (rji_length * rjk_length) - cos_theta * rji.x / (rji_length * rji_length))),
        static_cast<real_t>(f1 * (rjk.y / (rji_length * rjk_length) - cos_theta * rji.y / (rji_length * rji_length))),
        static_cast<real_t>(f1 * (rjk.z / (rji_length * rjk_length) - cos_theta * rji.z / (rji_length * rji_length)))};

    coord_t dk = {
        static_cast<real_t>(f1 * (rji.x / (rji_length * rjk_length) - cos_theta * rjk.x / (rjk_length * rjk_length))),
        static_cast<real_t>(f1 * (rji.y / (rji_length * rjk_length) - cos_theta * rjk.y / (rjk_length * rjk_length))),
        static_cast<real_t>(f1 * (rji.z / (rji_length * rjk_length) - cos_theta * rjk.z / (rjk_length * rjk_length)))};

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

real_t calc_angle_forces_host(int start, int end) {
    int N = end - start;
    if (N <= 0) return 0.0;
    using namespace CudaAngleForce;
    int blockSize = 256;
    int numBlocks = (N + blockSize - 1) / blockSize;

    auto &host_ctx = Context::instance();
    auto d_angles = host_ctx.angles->gpu_data_p;
    auto d_coords = host_ctx.coords->gpu_data_p;
    auto d_cangles = host_ctx.cangles->gpu_data_p;
    auto d_dvelocities = host_ctx.dvelocities->gpu_data_p;
    // todo: now have to do that, after moving all to CudaContext, can remove it
    // ctx.sync_all_to_device();

    real_t h_energy_sum = 0.0;
    cudaMemcpy(d_energy_sum, &h_energy_sum, sizeof(real_t), cudaMemcpyHostToDevice);

    // launch kernel
    calc_angle_forces_kernel<<<numBlocks, blockSize>>>(start, end, d_angles, d_coords, d_cangles, d_dvelocities, d_energy_sum);
    cudaDeviceSynchronize();

    // todo: Now have to do that, after moving all to CudaContext, can remove it
    // copy results back to host
    cudaMemcpy(&h_energy_sum, d_energy_sum, sizeof(real_t), cudaMemcpyDeviceToHost);
    return h_energy_sum;
}

void init_angle_force_kernel_data() {
    using namespace CudaAngleForce;
    if (!is_initialized) {
        check_cudaMalloc((void**)&d_energy_sum, sizeof(real_t));
        is_initialized = true;
    }
}

void cleanup_angle_force() {
    using namespace CudaAngleForce;
    if (is_initialized) {
        cudaFree(d_energy_sum);
        is_initialized = false;
    }
}
