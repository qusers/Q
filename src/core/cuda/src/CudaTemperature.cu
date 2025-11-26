#include "cuda/include/CudaTemperature.cuh"
#include "cuda/include/CudaUtility.cuh"
#include "iostream"
#include "utils.h"

namespace CudaTemperature {
bool is_initialized = false;
atype_t* d_atypes;
catype_t* d_catypes;
vel_t* d_velocities;
bool* d_excluded;
double* d_Temp_solute;
double* d_Tfree_solute;
double* d_Texcl_solute;
double* d_Temp_solvent;
double* d_Tfree_solvent;
double* d_Texcl_solvent;
}  // namespace CudaTemperature

__global__ void calc_temperature_kernel(int n_atoms, int n_atoms_solute, atype_t* atypes, catype_t* catypes, vel_t* velocities, bool* excluded, double boltz, double ekinmax,
                                        double* Temp_solute, double* Tfree_solute, double* Texcl_solute, double* Temp_solvent, double* Tfree_solvent, double* Texcl_solvent) {
    int idx = blockIdx.x * blockDim.x + threadIdx.x;
    if (idx >= n_atoms) return;
    double mass_i = catypes[atypes[idx].code - 1].m;
    double ener = .5 * mass_i * (pow(velocities[idx].x, 2) + pow(velocities[idx].y, 2) + pow(velocities[idx].z, 2));
    bool is_solute = (idx < n_atoms_solute);
    bool is_excluded = excluded[idx];

    if (is_solute) {
        atomicAdd(Temp_solute, ener);
        if (!is_excluded) {
            atomicAdd(Tfree_solute, ener);
        } else {
            atomicAdd(Texcl_solute, ener);
        }
    } else {
        atomicAdd(Temp_solvent, ener);
        if (!is_excluded) {
            atomicAdd(Tfree_solvent, ener);
        } else {
            atomicAdd(Texcl_solvent, ener);
        }
    }
    if (ener > ekinmax) {
        printf(">>> WARNING: hot atom %d: %f\n", idx, ener / Boltz / 3);
    }
}

void calc_temperature_host() {
    printf("Ndegf = %f, Ndegfree = %f, n_excluded = %d, Ndegfree_solvent = %f, Ndegfree_solute = %f\n", Ndegf, Ndegfree, n_excluded, Ndegfree_solvent, Ndegfree_solute);
    using namespace CudaTemperature;
    if (!is_initialized) {
        check_cudaMalloc((void**)&d_atypes, sizeof(atype_t) * n_atypes);
        check_cudaMalloc((void**)&d_catypes, sizeof(catype_t) * n_catypes);
        check_cudaMalloc((void**)&d_velocities, sizeof(vel_t) * n_atoms);
        check_cudaMalloc((void**)&d_excluded, sizeof(bool) * n_atoms);
        check_cudaMalloc((void**)&d_Temp_solute, sizeof(double));
        check_cudaMalloc((void**)&d_Tfree_solute, sizeof(double));
        check_cudaMalloc((void**)&d_Texcl_solute, sizeof(double));
        check_cudaMalloc((void**)&d_Temp_solvent, sizeof(double));
        check_cudaMalloc((void**)&d_Tfree_solvent, sizeof(double));
        check_cudaMalloc((void**)&d_Texcl_solvent, sizeof(double));

        cudaMemcpy(d_atypes, atypes, sizeof(atype_t) * n_atypes, cudaMemcpyHostToDevice);
        cudaMemcpy(d_catypes, catypes, sizeof(catype_t) * n_catypes, cudaMemcpyHostToDevice);
        cudaMemcpy(d_excluded, excluded, sizeof(bool) * n_atoms, cudaMemcpyHostToDevice);

        is_initialized = true;
    }
    cudaMemcpy(d_velocities, velocities, sizeof(vel_t) * n_atoms, cudaMemcpyHostToDevice);

    double h_Temp_solute = 0.0, h_Tfree_solute = 0.0, h_Texcl_solute = 0.0, h_Temp_solvent = 0.0, h_Tfree_solvent = 0.0, h_Texcl_solvent = 0.0;

    cudaMemcpy(d_Temp_solute, &h_Temp_solute, sizeof(double), cudaMemcpyHostToDevice);
    cudaMemcpy(d_Tfree_solute, &h_Tfree_solute, sizeof(double), cudaMemcpyHostToDevice);
    cudaMemcpy(d_Texcl_solute, &h_Texcl_solute, sizeof(double), cudaMemcpyHostToDevice);
    cudaMemcpy(d_Temp_solvent, &h_Temp_solvent, sizeof(double), cudaMemcpyHostToDevice);
    cudaMemcpy(d_Tfree_solvent, &h_Tfree_solvent, sizeof(double), cudaMemcpyHostToDevice);
    cudaMemcpy(d_Texcl_solvent, &h_Texcl_solvent, sizeof(double), cudaMemcpyHostToDevice);


    int blockSize = 256;
    int numBlocks = (n_atoms + blockSize - 1) / blockSize;

    double Ekinmax = 1000.0 * Ndegf * Boltz * md.temperature / 2.0 / n_atoms;
    calc_temperature_kernel<<<numBlocks, blockSize>>>(n_atoms, n_atoms_solute, d_atypes, d_catypes, d_velocities, d_excluded, Boltz, Ekinmax,
                                                      d_Temp_solute, d_Tfree_solute, d_Texcl_solute, d_Temp_solvent, d_Tfree_solvent, d_Texcl_solvent);

    cudaDeviceSynchronize();
    cudaMemcpy(&h_Temp_solute, d_Temp_solute, sizeof(double), cudaMemcpyDeviceToHost);
    cudaMemcpy(&h_Tfree_solute, d_Tfree_solute, sizeof(double), cudaMemcpyDeviceToHost);
    cudaMemcpy(&h_Texcl_solute, d_Texcl_solute, sizeof(double), cudaMemcpyDeviceToHost);
    cudaMemcpy(&h_Temp_solvent, d_Temp_solvent, sizeof(double), cudaMemcpyDeviceToHost);
    cudaMemcpy(&h_Tfree_solvent, d_Tfree_solvent, sizeof(double), cudaMemcpyDeviceToHost);
    cudaMemcpy(&h_Texcl_solvent, d_Texcl_solvent, sizeof(double), cudaMemcpyDeviceToHost);
    Tfree = h_Tfree_solute + h_Tfree_solvent;
    Temp = h_Temp_solute + h_Temp_solvent;

    E_total.Ukin = Temp;

    Temp = 2.0 * Temp / Boltz / Ndegf;
    Tfree = 2.0 * Tfree / Boltz / Ndegfree;

    if (separate_scaling) {
        if (h_Tfree_solvent != 0) Tscale_solvent = sqrt(1 + (dt / tau_T) * (md.temperature / h_Tfree_solvent - 1.0));
        if (h_Tfree_solute != 0) Tscale_solute = sqrt(1 + (dt / tau_T) * (md.temperature / h_Tfree_solute - 1.0));
    } else {
        if (Tfree != 0) Tscale_solvent = sqrt(1 + (dt / tau_T) * (md.temperature / Tfree - 1.0));
        Tscale_solute = Tscale_solvent;
    }
    printf("Tscale = %f, tau_T = %f, Temp = %f, Tfree = %f\n", Tscale_solvent, tau_T, Temp, Tfree);
}

void cleanup_temperature() {
    using namespace CudaTemperature;
    if (is_initialized) {
        cudaFree(d_atypes);
        cudaFree(d_catypes);
        cudaFree(d_velocities);
        cudaFree(d_excluded);
        cudaFree(d_Temp_solute);
        cudaFree(d_Tfree_solute);
        cudaFree(d_Texcl_solute);
        cudaFree(d_Temp_solvent);
        cudaFree(d_Tfree_solvent);
        cudaFree(d_Texcl_solvent);
        is_initialized = false;
    }   
}
