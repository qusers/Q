#include <iostream>

#include "cuda/include/cuda_restrwall_force.cuh"
#include "common/include/context.h"
#include "cuda_force_accumulation.cuh"

namespace CudaRestrwallForce {
bool is_initialized = false;
energy_accum_t* d_energies;
}  // namespace CudaRestrwallForce

__global__ void calc_restrwall_forces_kernel(
    restrwall_t* restrwalls,
    int n_restrwalls,
    coord_t* coords,
    energy_accum_t* energies,
    dvel_t* dvelocities,
    bool* heavy, topo_t topo) {
    int idx = blockIdx.x * blockDim.x + threadIdx.x;
    if (idx >= n_restrwalls) return;

    real_t k, b, db, ener, dv, fexp;
    coord_t dr;

    int ir = idx;
    k = restrwalls[ir].k;
    for (int i = restrwalls[ir].ai - 1; i < restrwalls[ir].aj - 1; i++) {
        if (heavy[i] || restrwalls[ir].ih) {
            dr.x = coords[i].x - topo.solvent_center.x;
            dr.y = coords[i].y - topo.solvent_center.y;
            dr.z = coords[i].z - topo.solvent_center.z;

            b = sqrt(dr.x * dr.x + dr.y * dr.y + dr.z * dr.z);
            db = b - restrwalls[ir].d;

            if (db > 0) {
                ener = .5 * k * db * db - restrwalls[ir].dMorse;
                dv = k * db / b;
            } else {
                fexp = exp(restrwalls[ir].aMorse * db);
                ener = restrwalls[ir].dMorse * (fexp * fexp - 2 * fexp);
                dv = -2 * restrwalls[ir].dMorse * restrwalls[ir].aMorse * (fexp - fexp * fexp) / b;
            }

            atomic_add_energy(energies, ener);

            atomic_add_force(&dvelocities[i].x, dv * dr.x);
            atomic_add_force(&dvelocities[i].y, dv * dr.y);
            atomic_add_force(&dvelocities[i].z, dv * dr.z);
        }
    }
}

void calc_restrwall_forces_host() {
    auto& host = Context::instance();
    if (host.n_restrwalls == 0) return;
    using namespace CudaRestrwallForce;
    auto d_restrwalls = host.restrwalls->gpu_data_p;
    auto d_coords = host.coords->gpu_data_p;
    auto d_dvelocities = host.dvelocities->gpu_data_p;
    auto d_heavy = host.heavy->gpu_data_p;
    cudaMemset(d_energies, 0, sizeof(energy_accum_t));

    int blockSize = 256;
    int numBlocks = (host.n_restrwalls + blockSize - 1) / blockSize;
    calc_restrwall_forces_kernel<<<numBlocks, blockSize>>>(
        d_restrwalls,
        host.n_restrwalls,
        d_coords,
        d_energies,
        d_dvelocities, d_heavy, host.topo);
    cudaDeviceSynchronize();
    energy_accum_t h_energy;
    cudaMemcpy(&h_energy, d_energies, sizeof(energy_accum_t), cudaMemcpyDeviceToHost);
    const real_t energy = energy_from_accum(h_energy);
    printf("Restrwall energy: %f\n", energy);
    host.E_restraint.Upres += energy;
}

void init_restrwall_force_kernel_data() {
    using namespace CudaRestrwallForce;
    if (!is_initialized) {
        check_cudaMalloc((void**)&d_energies, sizeof(energy_accum_t));
        is_initialized = true;
    }
}

void cleanup_restrwall_force() {
    using namespace CudaRestrwallForce;
    if (is_initialized) {
        cudaFree(d_energies);
        is_initialized = false;
    }
}
