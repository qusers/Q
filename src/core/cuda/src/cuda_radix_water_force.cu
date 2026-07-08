#include <stdexcept>

#include "temperature.h"
#include "common/include/constants.h"
#include "common/include/context.h"
#include "cuda/include/cuda_radix_water_force.cuh"
#include "cuda_force_accumulation.cuh"
namespace CudaRadixWaterForce {
}  // namespace CudaRadixWaterForce

__global__ void calc_radix_water_forces_kernel(
    coord_t* coords,
    const double* temperature_results,
    int n_atoms_solute,
    int n_atoms,
    topo_t topo,
    md_t md,
    double Dwmz,
    double awmz,
    dvel_t* dvelocities,
    energy_accum_t* energy) {
    int i = blockIdx.x * blockDim.x + threadIdx.x;
    i = n_atoms_solute + i * 3;  // Process only oxygen atoms of water molecules
    if (i >= n_atoms) return;

    const double shift = (md.radial_force != 0.0)
                             ? sqrt(Boltz * temperature_results[R_TFREE] / md.radial_force)
                             : 0.0;

    const double dx = coords[i].x - topo.solvent_center.x;
    const double dy = coords[i].y - topo.solvent_center.y;
    const double dz = coords[i].z - topo.solvent_center.z;
    const double b = sqrt(dx * dx + dy * dy + dz * dz);
    const double db = b - (topo.solvent_radius - shift);

    double ener;
    double dv;
    if (db > 0.0) {
        ener = 0.5 * md.radial_force * db * db - Dwmz;
        dv = md.radial_force * db / b;
    } else {
        if (b > 0.0) {
            const double fexp = exp(awmz * db);
            ener = Dwmz * (fexp * fexp - 2.0 * fexp);
            dv = -2.0 * Dwmz * awmz * (fexp - fexp * fexp) / b;
        } else {
            dv = 0.0;
            ener = 0.0;
        }
    }

    // Update energy and forces
    atomic_add_energy(&energy[E_RESTR_RADX], ener);
    atomic_add_force(&dvelocities[i].x, dv * dx);
    atomic_add_force(&dvelocities[i].y, dv * dy);
    atomic_add_force(&dvelocities[i].z, dv * dz);
}

void calc_radix_water_forces_host(const double* temperature_results) {
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

    calc_radix_water_forces_kernel<<<numBlocks, blockSize>>>(d_coords, temperature_results,
                                                             host.n_atoms_solute,
                                                             host.n_atoms,
                                                             host.topo,
                                                             host.md,
                                                             host.Dwmz,
                                                             host.awmz,
                                                             d_dvelocities,
                                                             host.energy.device());
}

void init_radix_water_force_kernel_data() {
}

void cleanup_radix_water_force() {
}
