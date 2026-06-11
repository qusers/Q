#include <stdexcept>

#include "common/include/constants.h"
#include "common/include/context.h"
#include "cuda/include/cuda_radix_water_force.cuh"
#include "cuda/include/cuda_utility.cuh"
#include "cuda_force_accumulation.cuh"
namespace CudaRadixWaterForce {
bool is_initialized = false;
energy_accum_t* d_energy;
}  // namespace CudaRadixWaterForce

__global__ void calc_radix_water_forces_kernel(
    coord_t* coords,
    real_t shift,
    int n_atoms_solute,
    int n_atoms,
    topo_t topo,
    md_t md,
    real_t Dwmz,
    real_t awmz,
    dvel_t* dvelocities,
    energy_accum_t* energy) {
    int i = blockIdx.x * blockDim.x + threadIdx.x;
    i = n_atoms_solute + i * 3;  // Process only oxygen atoms of water molecules
    if (i >= n_atoms) return;

    coord_t dr;

    dr.x = coords[i].x - topo.solvent_center.x;
    dr.y = coords[i].y - topo.solvent_center.y;
    dr.z = coords[i].z - topo.solvent_center.z;
    real_t b = sqrt(dr.x * dr.x + dr.y * dr.y + dr.z * dr.z);
    real_t db = b - (topo.solvent_radius - shift);

    real_t ener, dv;
    if (db > 0) {
        ener = 0.5 * md.radial_force * db * db - Dwmz;
        dv = md.radial_force * db / b;
    } else {
        if (b > 0.0) {
            real_t fexp = exp(awmz * db);
            ener = Dwmz * (fexp * fexp - 2 * fexp);
            dv = -2 * Dwmz * awmz * (fexp - fexp * fexp) / b;
        } else {
            dv = 0;
            ener = 0;
        }
    }

    // Update energy and forces
    atomic_add_energy(energy, ener);
    atomic_add_force(&dvelocities[i].x, dv * dr.x);
    atomic_add_force(&dvelocities[i].y, dv * dr.y);
    atomic_add_force(&dvelocities[i].z, dv * dr.z);
}

void calc_radix_water_forces_host() {
    auto& host = Context::instance();
    int water_atoms = host.n_atoms - host.n_atoms_solute;
    if (water_atoms == 0) {
        return;
    }
    using namespace CudaRadixWaterForce;
    int blockSize = 256;
    if (water_atoms % 3 != 0) {
        throw std::runtime_error("Number of water atoms is not a multiple of 3");
    }
    int oxygen_atoms = water_atoms / 3;
    int numBlocks = (oxygen_atoms + blockSize - 1) / blockSize;

    auto d_coords = host.coords->gpu_data_p;
    auto d_dvelocities = host.dvelocities->gpu_data_p;
    check_cuda(cudaMemset(d_energy, 0, sizeof(energy_accum_t)));

    real_t shift;
    if (host.md.radial_force != 0) {
        shift = sqrt(Boltz * host.Tfree / host.md.radial_force);
    } else {
        shift = 0;
    }

    energy_accum_t energy = 0;
    calc_radix_water_forces_kernel<<<numBlocks, blockSize>>>(d_coords,
                                                             shift,
                                                             host.n_atoms_solute,
                                                             host.n_atoms,
                                                             host.topo,
                                                             host.md,
                                                             host.Dwmz,
                                                             host.awmz,
                                                             d_dvelocities,
                                                             d_energy);
    check_cuda(cudaDeviceSynchronize());
    check_cuda(cudaMemcpy(&energy, d_energy, sizeof(energy_accum_t), cudaMemcpyDeviceToHost));
    host.E_restraint.Uradx += energy_from_accum(energy);
}

void init_radix_water_force_kernel_data() {
    using namespace CudaRadixWaterForce;
    if (!is_initialized) {
        check_cudaMalloc((void**)&d_energy, sizeof(energy_accum_t));
        is_initialized = true;
    }
}

void cleanup_radix_water_force() {
    using namespace CudaRadixWaterForce;
    if (is_initialized) {
        cudaFree(d_energy);
        is_initialized = false;
    }
}
