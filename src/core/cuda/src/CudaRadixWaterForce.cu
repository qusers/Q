#include <stdexcept>

#include "cuda/include/CudaRadixWaterForce.cuh"
#include "cuda/include/CudaUtility.cuh"
#include "utils.h"

namespace CudaRadixWaterForce {
bool is_initialized = false;
coord_t* d_coords = nullptr;
dvel_t* d_dvelocities = nullptr;
double* d_energy = nullptr;
}  // namespace CudaRadixWaterForce
__global__ void calc_radix_water_forces_kernel(
    coord_t* coords, 
    double shift, 
    int n_atoms_solute,
    int n_atoms,
    topo_t topo,
    md_t md,
    double Dwmz,
    double awmz,
    dvel_t* dvelocities, 
    double* energy) {
    int i = blockIdx.x * blockDim.x + threadIdx.x;
    i = n_atoms_solute + i * 3;  // Process only oxygen atoms of water molecules
    if (i >= n_atoms) return;

    coord_t dr;

    dr.x = coords[i].x - topo.solvent_center.x;
    dr.y = coords[i].y - topo.solvent_center.y;
    dr.z = coords[i].z - topo.solvent_center.z;
    double b = sqrt(pow(dr.x, 2) + pow(dr.y, 2) + pow(dr.z, 2));
    double db = b - (topo.solvent_radius - shift);

    double ener, dv;
    if (db > 0) {
        ener = 0.5 * md.radial_force * pow(db, 2) - Dwmz;
        dv = md.radial_force * db / b;
    } else {
        if (b > 0.0) {
            double fexp = exp(awmz * db);
            ener = Dwmz * (pow(fexp, 2) - 2 * fexp);
            dv = -2 * Dwmz * awmz * (fexp - pow(fexp, 2)) / b;
        } else {
            dv = 0;
            ener = 0;
        }
    }

    // Update energy and forces
    atomicAdd(energy, ener);
    atomicAdd(&dvelocities[i].x, dv * dr.x);
    atomicAdd(&dvelocities[i].y, dv * dr.y);
    atomicAdd(&dvelocities[i].z, dv * dr.z);
}

void calc_radix_water_forces_host() {
    int water_atoms = n_atoms - n_atoms_solute;
    if (water_atoms == 0) {
        return;
    }
    int blockSize = 256;
    if (water_atoms % 3 != 0) {
        throw std::runtime_error("Number of water atoms is not a multiple of 3");
    }
    int oxygen_atoms = water_atoms / 3;
    int numBlocks = (oxygen_atoms + blockSize - 1) / blockSize;

    using namespace CudaRadixWaterForce;
    if (!is_initialized) {
        check_cudaMalloc((void**)&d_coords, sizeof(coord_t) * n_atoms);
        check_cudaMalloc((void**)&d_dvelocities, sizeof(dvel_t) * n_atoms);
        check_cudaMalloc((void**)&d_energy, sizeof(double));
        is_initialized = true;
    }

    cudaMemcpy(d_coords, coords, sizeof(coord_t) * n_atoms, cudaMemcpyHostToDevice);
    cudaMemcpy(d_dvelocities, dvelocities, sizeof(dvel_t) * n_atoms, cudaMemcpyHostToDevice);
    double energy = 0.0;
    cudaMemcpy(d_energy, &energy, sizeof(double), cudaMemcpyHostToDevice);

    double shift;
    if (md.radial_force != 0) {
        shift = sqrt(Boltz * Tfree / md.radial_force);
    } else {
        shift = 0;
    }
    calc_radix_water_forces_kernel<<<numBlocks, blockSize>>>(d_coords, 
        shift, 
        n_atoms_solute, 
        n_atoms, 
        topo, 
        md, 
        Dwmz, 
        awmz, 
        d_dvelocities, 
        d_energy);
    cudaDeviceSynchronize();
    cudaMemcpy(dvelocities, d_dvelocities, sizeof(dvel_t) * n_atoms, cudaMemcpyDeviceToHost);
    cudaMemcpy(&energy, d_energy, sizeof(double), cudaMemcpyDeviceToHost);
    E_restraint.Uradx += energy;
}

void cleanup_radix_water_force() {
    using namespace CudaRadixWaterForce;
    if (is_initialized) {
        cudaFree(d_coords);
        cudaFree(d_dvelocities);
        cudaFree(d_energy);
        is_initialized = false;
    }
}
