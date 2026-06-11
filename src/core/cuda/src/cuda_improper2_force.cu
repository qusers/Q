#include "cuda/include/cuda_improper2_force.cuh"
#include "cuda/include/cuda_utility.cuh"
#include "context.h"
#include "cuda_force_accumulation.cuh"

namespace CudaImproper2Force {
bool is_initialized = false;
energy_accum_t* d_energy_sum;
}  // namespace CudaImproper2Force

__global__ void calc_improper2_forces_kernel(int start, int end, improper_t* impropers, cimproper_t* cimpropers, coord_t* coords, dvel_t* dvelocities, energy_accum_t* energy_sum) {
    int i = blockIdx.x * blockDim.x + threadIdx.x + start;
    if (i >= end) return;

    int aii, aji, aki, ali;

    coord_t ai, aj, ak, al;

    improper_t imp;
    cimproper_t cimp;

    imp = impropers[i];
    cimp = cimpropers[imp.code - 1];

    aii = imp.ai - 1;
    aji = imp.aj - 1;
    aki = imp.ak - 1;
    ali = imp.al - 1;

    ai = coords[aii];
    aj = coords[aji];
    ak = coords[aki];
    al = coords[ali];

    const double rji_x = static_cast<double>(ai.x) - static_cast<double>(aj.x);
    const double rji_y = static_cast<double>(ai.y) - static_cast<double>(aj.y);
    const double rji_z = static_cast<double>(ai.z) - static_cast<double>(aj.z);
    const double rjk_x = static_cast<double>(ak.x) - static_cast<double>(aj.x);
    const double rjk_y = static_cast<double>(ak.y) - static_cast<double>(aj.y);
    const double rjk_z = static_cast<double>(ak.z) - static_cast<double>(aj.z);
    const double rkl_x = static_cast<double>(al.x) - static_cast<double>(ak.x);
    const double rkl_y = static_cast<double>(al.y) - static_cast<double>(ak.y);
    const double rkl_z = static_cast<double>(al.z) - static_cast<double>(ak.z);

    const double rnj_x = rji_y * rjk_z - rji_z * rjk_y;
    const double rnj_y = rji_z * rjk_x - rji_x * rjk_z;
    const double rnj_z = rji_x * rjk_y - rji_y * rjk_x;
    const double rnk_x = -rjk_y * rkl_z + rjk_z * rkl_y;
    const double rnk_y = -rjk_z * rkl_x + rjk_x * rkl_z;
    const double rnk_z = -rjk_x * rkl_y + rjk_y * rkl_x;

    const double bj2inv = 1.0 / (rnj_x * rnj_x + rnj_y * rnj_y + rnj_z * rnj_z);
    const double bk2inv = 1.0 / (rnk_x * rnk_x + rnk_y * rnk_y + rnk_z * rnk_z);
    const double bjinv = sqrt(bj2inv);
    const double bkinv = sqrt(bk2inv);
    const double bjkinv = bjinv * bkinv;

    double cos_phi = (rnj_x * rnk_x + rnj_y * rnk_y + rnj_z * rnk_z) * bjkinv;
    // printf("cos_phi = %f\n", cos_phi);
    if (cos_phi > 1) {
        cos_phi = 1;
    }
    if (cos_phi < -1) {
        cos_phi = -1;
    }
    double phi = acos(cos_phi);
    if (rjk_x * (rnj_y * rnk_z - rnj_z * rnk_y) + rjk_y * (rnj_z * rnk_x - rnj_x * rnk_z) + rjk_z * (rnj_x * rnk_y - rnj_y * rnk_x) < 0) {
        phi = -phi;
    }

    // Energy
    const double arg = 2.0 * phi - to_radians_device(cimp.phi0);
    const double ener = cimp.k * (1.0 + cos(arg));
    const double dv = -2.0 * cimp.k * sin(arg);

    // Forces
    double f1 = sin(phi);
    const double sin_epsilon = static_cast<double>(k_singular_sin_epsilon);
    if (fabs(f1) < sin_epsilon) f1 = copysign(sin_epsilon, f1);
    f1 = -1.0 / f1;
    // printf("f1 = %f phi = %f cos_phi = %f\n", f1, phi, cos_phi);

    const double di_x = f1 * (rnk_x * bjkinv - cos_phi * rnj_x * bj2inv);
    const double di_y = f1 * (rnk_y * bjkinv - cos_phi * rnj_y * bj2inv);
    const double di_z = f1 * (rnk_z * bjkinv - cos_phi * rnj_z * bj2inv);
    const double dl_x = f1 * (rnj_x * bjkinv - cos_phi * rnk_x * bk2inv);
    const double dl_y = f1 * (rnj_y * bjkinv - cos_phi * rnk_y * bk2inv);
    const double dl_z = f1 * (rnj_z * bjkinv - cos_phi * rnk_z * bk2inv);

    const double rki_x = rji_x - rjk_x;
    const double rki_y = rji_y - rjk_y;
    const double rki_z = rji_z - rjk_z;
    const double rlj_x = -rjk_x - rkl_x;
    const double rlj_y = -rjk_y - rkl_y;
    const double rlj_z = -rjk_z - rkl_z;

    const double dpi_x = rjk_y * di_z - rjk_z * di_y;
    const double dpi_y = rjk_z * di_x - rjk_x * di_z;
    const double dpi_z = rjk_x * di_y - rjk_y * di_x;
    const double dpj_x = rki_y * di_z - rki_z * di_y + rkl_y * dl_z - rkl_z * dl_y;
    const double dpj_y = rki_z * di_x - rki_x * di_z + rkl_z * dl_x - rkl_x * dl_z;
    const double dpj_z = rki_x * di_y - rki_y * di_x + rkl_x * dl_y - rkl_y * dl_x;
    const double dpk_x = rlj_y * dl_z - rlj_z * dl_y - rji_y * di_z + rji_z * di_y;
    const double dpk_y = rlj_z * dl_x - rlj_x * dl_z - rji_z * di_x + rji_x * di_z;
    const double dpk_z = rlj_x * dl_y - rlj_y * dl_x - rji_x * di_y + rji_y * di_x;
    const double dpl_x = rjk_y * dl_z - rjk_z * dl_y;
    const double dpl_y = rjk_z * dl_x - rjk_x * dl_z;
    const double dpl_z = rjk_x * dl_y - rjk_y * dl_x;

    // Update energy and forces
    atomic_add_energy(energy_sum, ener);

    atomic_add_force(&dvelocities[aii].x, dv * dpi_x);
    atomic_add_force(&dvelocities[aii].y, dv * dpi_y);
    atomic_add_force(&dvelocities[aii].z, dv * dpi_z);
    atomic_add_force(&dvelocities[aji].x, dv * dpj_x);
    atomic_add_force(&dvelocities[aji].y, dv * dpj_y);
    atomic_add_force(&dvelocities[aji].z, dv * dpj_z);
    atomic_add_force(&dvelocities[aki].x, dv * dpk_x);
    atomic_add_force(&dvelocities[aki].y, dv * dpk_y);
    atomic_add_force(&dvelocities[aki].z, dv * dpk_z);
    atomic_add_force(&dvelocities[ali].x, dv * dpl_x);
    atomic_add_force(&dvelocities[ali].y, dv * dpl_y);
    atomic_add_force(&dvelocities[ali].z, dv * dpl_z);
}

real_t calc_improper2_forces_host(int start, int end) {
    int N = end - start;
    if (N <= 0) return 0.0;
    using namespace CudaImproper2Force;
    int blockSize = 256;
    int numBlocks = (N + blockSize - 1) / blockSize;

    energy_accum_t energy = 0;
    cudaMemcpy(d_energy_sum, &energy, sizeof(energy_accum_t), cudaMemcpyHostToDevice);

    auto& host_ctx = Context::instance();
    coord_t* d_coords = host_ctx.coords->gpu_data_p;
    dvel_t* d_dvelocities = host_ctx.dvelocities->gpu_data_p;
    improper_t* d_impropers = host_ctx.impropers->gpu_data_p;
    cimproper_t* d_cimpropers = host_ctx.cimpropers->gpu_data_p;

    calc_improper2_forces_kernel<<<numBlocks, blockSize>>>(start, end, d_impropers, d_cimpropers, d_coords, d_dvelocities, d_energy_sum);
    cudaDeviceSynchronize();
    cudaMemcpy(&energy, d_energy_sum, sizeof(energy_accum_t), cudaMemcpyDeviceToHost);
    return energy_from_accum(energy);
}

void init_improper2_force_kernel_data() {
    using namespace CudaImproper2Force;
    if (!is_initialized) {
        check_cudaMalloc((void**)&d_energy_sum, sizeof(energy_accum_t));
        is_initialized = true;
    }
}

void cleanup_improper2_force() {
    using namespace CudaImproper2Force;
    if (is_initialized) {
        cudaFree(d_energy_sum);
        is_initialized = false;
    }
}
