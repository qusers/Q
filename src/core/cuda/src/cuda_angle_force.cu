#include "cuda/include/cuda_angle_force.cuh"
#include "cuda/include/cuda_utility.cuh"
#include "context.h"
#include "cuda_force_accumulation.cuh"

namespace CudaAngleForce {
bool is_initialized = false;
energy_accum_t* d_energy_sum;
}  // namespace CudaAngleForce

__global__ void calc_angle_forces_kernel(int start, int end, angle_t* angles, coord_t* coords, cangle_t* cangles, dvel_t* dvelocities, energy_accum_t* energy_sum) {
    int idx = blockIdx.x * blockDim.x + threadIdx.x + start;
    if (idx >= end) return;

    int i = angles[idx].ai - 1;
    int j = angles[idx].aj - 1;
    int k = angles[idx].ak - 1;

    coord_t ri = coords[i];
    coord_t rj = coords[j];
    coord_t rk = coords[k];

    cangle_t cang = cangles[angles[idx].code - 1];

    const double rji_x = ri.x - rj.x;
    const double rji_y = ri.y - rj.y;
    const double rji_z = ri.z - rj.z;

    const double rjk_x = rk.x - rj.x;
    const double rjk_y = rk.y - rj.y;
    const double rjk_z = rk.z - rj.z;

    const double rji2inv = 1.0 / (rji_x * rji_x + rji_y * rji_y + rji_z * rji_z);
    const double rjk2inv = 1.0 / (rjk_x * rjk_x + rjk_y * rjk_y + rjk_z * rjk_z);
    const double rjiinv = sqrt(rji2inv);
    const double rjkinv = sqrt(rjk2inv);

    double cos_theta = (rji_x * rjk_x + rji_y * rjk_y + rji_z * rjk_z) * rjiinv * rjkinv;

    cos_theta = cos_theta > 1.0 ? 1.0 : cos_theta;
    cos_theta = cos_theta < -1.0 ? -1.0 : cos_theta;
    const double theta = acos(cos_theta);

    const double dtheta = theta - to_radians_device(cang.th0);
    const double energy = 0.5 * cang.kth * dtheta * dtheta;

    // calculate force magnitude
    const double dv = cang.kth * dtheta;

    double f1 = sin(theta);
    const double sin_epsilon = k_singular_sin_epsilon;
    if (fabs(f1) < sin_epsilon) {
        f1 = -1.0 / sin_epsilon;
    } else {
        f1 = -1.0 / f1;
    }

    atomic_add_energy(energy_sum, energy);

    const double di_x = f1 * (rjk_x * rjiinv * rjkinv - cos_theta * rji_x * rji2inv);
    const double di_y = f1 * (rjk_y * rjiinv * rjkinv - cos_theta * rji_y * rji2inv);
    const double di_z = f1 * (rjk_z * rjiinv * rjkinv - cos_theta * rji_z * rji2inv);

    const double dk_x = f1 * (rji_x * rjiinv * rjkinv - cos_theta * rjk_x * rjk2inv);
    const double dk_y = f1 * (rji_y * rjiinv * rjkinv - cos_theta * rjk_y * rjk2inv);
    const double dk_z = f1 * (rji_z * rjiinv * rjkinv - cos_theta * rjk_z * rjk2inv);

    atomic_add_force(&dvelocities[i].x, dv * di_x);
    atomic_add_force(&dvelocities[i].y, dv * di_y);
    atomic_add_force(&dvelocities[i].z, dv * di_z);

    atomic_add_force(&dvelocities[k].x, dv * dk_x);
    atomic_add_force(&dvelocities[k].y, dv * dk_y);
    atomic_add_force(&dvelocities[k].z, dv * dk_z);

    atomic_add_force(&dvelocities[j].x, -dv * (di_x + dk_x));
    atomic_add_force(&dvelocities[j].y, -dv * (di_y + dk_y));
    atomic_add_force(&dvelocities[j].z, -dv * (di_z + dk_z));
}

double calc_angle_forces_host(int start, int end) {
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

    energy_accum_t h_energy_sum = 0;
    cudaMemcpy(d_energy_sum, &h_energy_sum, sizeof(energy_accum_t), cudaMemcpyHostToDevice);

    // launch kernel
    calc_angle_forces_kernel<<<numBlocks, blockSize>>>(start, end, d_angles, d_coords, d_cangles, d_dvelocities, d_energy_sum);
    cudaDeviceSynchronize();

    // todo: Now have to do that, after moving all to CudaContext, can remove it
    // copy results back to host
    cudaMemcpy(&h_energy_sum, d_energy_sum, sizeof(energy_accum_t), cudaMemcpyDeviceToHost);
    return energy_from_accum(h_energy_sum);
}

void init_angle_force_kernel_data() {
    using namespace CudaAngleForce;
    if (!is_initialized) {
        check_cudaMalloc((void**)&d_energy_sum, sizeof(energy_accum_t));
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
