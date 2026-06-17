#include "cuda/include/cuda_temperature.cuh"
#include "cuda/include/cuda_utility.cuh"
#include "common/include/constants.h"
#include "common/include/context.h"
#include "cuda_force_accumulation.cuh"
#include "iostream"

namespace CudaTemperature {
bool is_initialized = false;
energy_accum_t* d_Temp_solute;
energy_accum_t* d_Tfree_solute;
energy_accum_t* d_Texcl_solute;
energy_accum_t* d_Temp_solvent;
energy_accum_t* d_Tfree_solvent;
energy_accum_t* d_Texcl_solvent;
}  // namespace CudaTemperature

__global__ void calc_temperature_kernel(int n_atoms, int n_atoms_solute, atype_t* atypes, catype_t* catypes, vel_t* velocities, bool* excluded, double boltz, double ekinmax,
                                        energy_accum_t* Temp_solute, energy_accum_t* Tfree_solute, energy_accum_t* Texcl_solute, energy_accum_t* Temp_solvent, energy_accum_t* Tfree_solvent, energy_accum_t* Texcl_solvent) {
    int idx = blockIdx.x * blockDim.x + threadIdx.x;
    if (idx >= n_atoms) return;
    double mass_i = catypes[atypes[idx].code - 1].m;
    const double vx = velocities[idx].x;
    const double vy = velocities[idx].y;
    const double vz = velocities[idx].z;
    double ener = .5 * mass_i * (vx * vx + vy * vy + vz * vz);
    bool is_solute = (idx < n_atoms_solute);
    bool is_excluded = excluded[idx];

    if (is_solute) {
        atomic_add_energy(Temp_solute, ener);
        if (!is_excluded) {
            atomic_add_energy(Tfree_solute, ener);
        } else {
            atomic_add_energy(Texcl_solute, ener);
        }
    } else {
        atomic_add_energy(Temp_solvent, ener);
        if (!is_excluded) {
            atomic_add_energy(Tfree_solvent, ener);
        } else {
            atomic_add_energy(Texcl_solvent, ener);
        }
    }
    if (ener > ekinmax) {
        printf(">>> WARNING: hot atom %d: %f\n", idx, ener / boltz / 3);
    }
}

void calc_temperature_host() {
    auto& host = Context::instance();
    using namespace CudaTemperature;
    energy_accum_t h_Temp_solute = 0, h_Tfree_solute = 0, h_Texcl_solute = 0, h_Temp_solvent = 0, h_Tfree_solvent = 0, h_Texcl_solvent = 0;

    cudaMemcpy(d_Temp_solute, &h_Temp_solute, sizeof(energy_accum_t), cudaMemcpyHostToDevice);
    cudaMemcpy(d_Tfree_solute, &h_Tfree_solute, sizeof(energy_accum_t), cudaMemcpyHostToDevice);
    cudaMemcpy(d_Texcl_solute, &h_Texcl_solute, sizeof(energy_accum_t), cudaMemcpyHostToDevice);
    cudaMemcpy(d_Temp_solvent, &h_Temp_solvent, sizeof(energy_accum_t), cudaMemcpyHostToDevice);
    cudaMemcpy(d_Tfree_solvent, &h_Tfree_solvent, sizeof(energy_accum_t), cudaMemcpyHostToDevice);
    cudaMemcpy(d_Texcl_solvent, &h_Texcl_solvent, sizeof(energy_accum_t), cudaMemcpyHostToDevice);

    atype_t* d_atypes = host.atypes->gpu_data_p;
    catype_t* d_catypes = host.catypes->gpu_data_p;
    vel_t* d_velocities = host.velocities->gpu_data_p;
    bool* d_excluded = host.excluded->gpu_data_p;

    int blockSize = 256;
    int numBlocks = (host.n_atoms + blockSize - 1) / blockSize;

    double Ekinmax = 1000.0 * host.Ndegf * Boltz * host.md.temperature / 2.0 / host.n_atoms;
    calc_temperature_kernel<<<numBlocks, blockSize>>>(host.n_atoms, host.n_atoms_solute, d_atypes, d_catypes, d_velocities, d_excluded, Boltz, Ekinmax,
                                                      d_Temp_solute, d_Tfree_solute, d_Texcl_solute, d_Temp_solvent, d_Tfree_solvent, d_Texcl_solvent);

    cudaDeviceSynchronize();
    cudaMemcpy(&h_Temp_solute, d_Temp_solute, sizeof(energy_accum_t), cudaMemcpyDeviceToHost);
    cudaMemcpy(&h_Tfree_solute, d_Tfree_solute, sizeof(energy_accum_t), cudaMemcpyDeviceToHost);
    cudaMemcpy(&h_Texcl_solute, d_Texcl_solute, sizeof(energy_accum_t), cudaMemcpyDeviceToHost);
    cudaMemcpy(&h_Temp_solvent, d_Temp_solvent, sizeof(energy_accum_t), cudaMemcpyDeviceToHost);
    cudaMemcpy(&h_Tfree_solvent, d_Tfree_solvent, sizeof(energy_accum_t), cudaMemcpyDeviceToHost);
    cudaMemcpy(&h_Texcl_solvent, d_Texcl_solvent, sizeof(energy_accum_t), cudaMemcpyDeviceToHost);
    double Tfree_solute_value = energy_from_accum(h_Tfree_solute);
    double Tfree_solvent_value = energy_from_accum(h_Tfree_solvent);
    host.Tfree = Tfree_solute_value + Tfree_solvent_value;
    host.Temp = energy_from_accum(h_Temp_solute) + energy_from_accum(h_Temp_solvent);

    host.E_total.Ukin = host.Temp;

    host.Temp = 2.0 * host.Temp / Boltz / host.Ndegf;
    host.Tfree = 2.0 * host.Tfree / Boltz / host.Ndegfree;

    if (host.separate_scaling) {
        Tfree_solvent_value = 2.0 * Tfree_solvent_value / Boltz / host.Ndegfree_solvent;
        Tfree_solute_value = 2.0 * Tfree_solute_value / Boltz / host.Ndegfree_solute;
        if (Tfree_solvent_value != 0) host.Tscale_solvent = sqrt(1 + (host.dt / host.tau_T) * (host.md.temperature / Tfree_solvent_value - 1.0));
        if (Tfree_solute_value != 0) host.Tscale_solute = sqrt(1 + (host.dt / host.tau_T) * (host.md.temperature / Tfree_solute_value - 1.0));
    } else {
        if (host.Tfree != 0) host.Tscale_solvent = sqrt(1 + (host.dt / host.tau_T) * (host.md.temperature / host.Tfree - 1.0));
        host.Tscale_solute = host.Tscale_solvent;
    }
    printf("Tscale = %f, tau_T = %f, Temp = %f, Tfree = %f\n", host.Tscale_solvent, host.tau_T, host.Temp, host.Tfree);
}

void init_temperature_kernel_data() {
    using namespace CudaTemperature;
    if (!is_initialized) {
        check_cudaMalloc((void**)&d_Temp_solute, sizeof(energy_accum_t));
        check_cudaMalloc((void**)&d_Tfree_solute, sizeof(energy_accum_t));
        check_cudaMalloc((void**)&d_Texcl_solute, sizeof(energy_accum_t));
        check_cudaMalloc((void**)&d_Temp_solvent, sizeof(energy_accum_t));
        check_cudaMalloc((void**)&d_Tfree_solvent, sizeof(energy_accum_t));
        check_cudaMalloc((void**)&d_Texcl_solvent, sizeof(energy_accum_t));
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
