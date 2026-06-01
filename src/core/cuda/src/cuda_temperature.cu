#include "cuda/include/cuda_temperature.cuh"
#include "cuda/include/cuda_force_accum.cuh"
#include "cuda/include/cuda_utility.cuh"
#include "common/include/constants.h"
#include "common/include/context.h"
#include "iostream"

namespace CudaTemperature {
bool is_initialized = false;
force_fixed_storage_t* d_Temp_solute;
force_fixed_storage_t* d_Tfree_solute;
force_fixed_storage_t* d_Texcl_solute;
force_fixed_storage_t* d_Temp_solvent;
force_fixed_storage_t* d_Tfree_solvent;
force_fixed_storage_t* d_Texcl_solvent;
}  // namespace CudaTemperature

__global__ void calc_temperature_kernel(int n_atoms, int n_atoms_solute, atype_t* atypes, catype_t* catypes, vel_t* velocities, bool* excluded, real_t boltz, real_t ekinmax,
                                        force_fixed_storage_t* Temp_solute, force_fixed_storage_t* Tfree_solute, force_fixed_storage_t* Texcl_solute,
                                        force_fixed_storage_t* Temp_solvent, force_fixed_storage_t* Tfree_solvent, force_fixed_storage_t* Texcl_solvent) {
    int idx = blockIdx.x * blockDim.x + threadIdx.x;
    if (idx >= n_atoms) return;
    real_t mass_i = catypes[atypes[idx].code - 1].m;
    const real_t vx = velocities[idx].x;
    const real_t vy = velocities[idx].y;
    const real_t vz = velocities[idx].z;
    real_t ener = .5 * mass_i * (vx * vx + vy * vy + vz * vz);
    bool is_solute = (idx < n_atoms_solute);
    bool is_excluded = excluded[idx];

    if (is_solute) {
        atomic_add_force_component(Temp_solute, ener);
        if (!is_excluded) {
            atomic_add_force_component(Tfree_solute, ener);
        } else {
            atomic_add_force_component(Texcl_solute, ener);
        }
    } else {
        atomic_add_force_component(Temp_solvent, ener);
        if (!is_excluded) {
            atomic_add_force_component(Tfree_solvent, ener);
        } else {
            atomic_add_force_component(Texcl_solvent, ener);
        }
    }
    if (ener > ekinmax) {
        printf(">>> WARNING: hot atom %d: %f\n", idx, ener / boltz / 3);
    }
}

void calc_temperature_host() {
    auto& host = Context::instance();
    using namespace CudaTemperature;
    real_t h_Temp_solute = 0.0, h_Tfree_solute = 0.0, h_Temp_solvent = 0.0, h_Tfree_solvent = 0.0;

    cudaMemset(d_Temp_solute, 0, sizeof(force_fixed_storage_t));
    cudaMemset(d_Tfree_solute, 0, sizeof(force_fixed_storage_t));
    cudaMemset(d_Texcl_solute, 0, sizeof(force_fixed_storage_t));
    cudaMemset(d_Temp_solvent, 0, sizeof(force_fixed_storage_t));
    cudaMemset(d_Tfree_solvent, 0, sizeof(force_fixed_storage_t));
    cudaMemset(d_Texcl_solvent, 0, sizeof(force_fixed_storage_t));

    atype_t* d_atypes = host.atypes->gpu_data_p;
    catype_t* d_catypes = host.catypes->gpu_data_p;
    vel_t* d_velocities = host.velocities->gpu_data_p;
    bool* d_excluded = host.excluded->gpu_data_p;

    int blockSize = 256;
    int numBlocks = (host.n_atoms + blockSize - 1) / blockSize;

    real_t Ekinmax = 1000.0 * host.Ndegf * Boltz * host.md.temperature / 2.0 / host.n_atoms;
    calc_temperature_kernel<<<numBlocks, blockSize>>>(host.n_atoms, host.n_atoms_solute, d_atypes, d_catypes, d_velocities, d_excluded, Boltz, Ekinmax,
                                                      d_Temp_solute, d_Tfree_solute, d_Texcl_solute, d_Temp_solvent, d_Tfree_solvent, d_Texcl_solvent);

    cudaDeviceSynchronize();
    force_fixed_storage_t temp_solute_fixed = 0;
    force_fixed_storage_t tfree_solute_fixed = 0;
    force_fixed_storage_t temp_solvent_fixed = 0;
    force_fixed_storage_t tfree_solvent_fixed = 0;
    cudaMemcpy(&temp_solute_fixed, d_Temp_solute, sizeof(force_fixed_storage_t), cudaMemcpyDeviceToHost);
    cudaMemcpy(&tfree_solute_fixed, d_Tfree_solute, sizeof(force_fixed_storage_t), cudaMemcpyDeviceToHost);
    cudaMemcpy(&temp_solvent_fixed, d_Temp_solvent, sizeof(force_fixed_storage_t), cudaMemcpyDeviceToHost);
    cudaMemcpy(&tfree_solvent_fixed, d_Tfree_solvent, sizeof(force_fixed_storage_t), cudaMemcpyDeviceToHost);
    h_Temp_solute = fixed_to_force(force_fixed_from_storage(temp_solute_fixed));
    h_Tfree_solute = fixed_to_force(force_fixed_from_storage(tfree_solute_fixed));
    h_Temp_solvent = fixed_to_force(force_fixed_from_storage(temp_solvent_fixed));
    h_Tfree_solvent = fixed_to_force(force_fixed_from_storage(tfree_solvent_fixed));
    host.Tfree = h_Tfree_solute + h_Tfree_solvent;
    host.Temp = h_Temp_solute + h_Temp_solvent;

    host.E_total.Ukin = host.Temp;

    host.Temp = 2.0 * host.Temp / Boltz / host.Ndegf;
    host.Tfree = 2.0 * host.Tfree / Boltz / host.Ndegfree;

    if (host.separate_scaling) {
        h_Tfree_solvent = 2.0 * h_Tfree_solvent / Boltz / host.Ndegfree_solvent;
        h_Tfree_solute = 2.0 * h_Tfree_solute / Boltz / host.Ndegfree_solute;
        if (h_Tfree_solvent != 0) host.Tscale_solvent = sqrt(1 + (host.dt / host.tau_T) * (host.md.temperature / h_Tfree_solvent - 1.0));
        if (h_Tfree_solute != 0) host.Tscale_solute = sqrt(1 + (host.dt / host.tau_T) * (host.md.temperature / h_Tfree_solute - 1.0));
    } else {
        if (host.Tfree != 0) host.Tscale_solvent = sqrt(1 + (host.dt / host.tau_T) * (host.md.temperature / host.Tfree - 1.0));
        host.Tscale_solute = host.Tscale_solvent;
    }
    printf("Tscale = %f, tau_T = %f, Temp = %f, Tfree = %f\n", host.Tscale_solvent, host.tau_T, host.Temp, host.Tfree);
}

void init_temperature_kernel_data() {
    using namespace CudaTemperature;
    if (!is_initialized) {
        check_cudaMalloc((void**)&d_Temp_solute, sizeof(force_fixed_storage_t));
        check_cudaMalloc((void**)&d_Tfree_solute, sizeof(force_fixed_storage_t));
        check_cudaMalloc((void**)&d_Texcl_solute, sizeof(force_fixed_storage_t));
        check_cudaMalloc((void**)&d_Temp_solvent, sizeof(force_fixed_storage_t));
        check_cudaMalloc((void**)&d_Tfree_solvent, sizeof(force_fixed_storage_t));
        check_cudaMalloc((void**)&d_Texcl_solvent, sizeof(force_fixed_storage_t));
        is_initialized = true;
    }
}

void cleanup_temperature() {
    using namespace CudaTemperature;
    if (is_initialized) {
        cudaFree(d_Temp_solute);
        cudaFree(d_Tfree_solute);
        cudaFree(d_Texcl_solute);
        cudaFree(d_Temp_solvent);
        cudaFree(d_Tfree_solvent);
        cudaFree(d_Texcl_solvent);
        is_initialized = false;
    }
}
