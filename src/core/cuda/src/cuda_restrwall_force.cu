#include <iostream>

#include "cuda/include/cuda_restrwall_force.cuh"
#include "common/include/context.h"
#include "cuda_force_accumulation.cuh"

namespace CudaRestrwallForce {
}  // namespace CudaRestrwallForce

__global__ void calc_restrwall_forces_kernel(
    restrwall_t* restrwalls,
    int n_restrwalls,
    coord_t* coords,
    energy_accum_t* energy,
    dvel_t* dvelocities,
    bool* heavy, topo_t topo) {
    int idx = blockIdx.x * blockDim.x + threadIdx.x;
    if (idx >= n_restrwalls) return;

    int ir = idx;
    const double k = restrwalls[ir].k;
    for (int i = restrwalls[ir].ai - 1; i < restrwalls[ir].aj - 1; i++) {
        if (heavy[i] || restrwalls[ir].ih) {
            const double dx = coords[i].x - topo.solvent_center.x;
            const double dy = coords[i].y - topo.solvent_center.y;
            const double dz = coords[i].z - topo.solvent_center.z;

            const double b = sqrt(dx * dx + dy * dy + dz * dz);
            const double db = b - restrwalls[ir].d;

            double ener;
            double dv;
            if (db > 0.0) {
                ener = 0.5 * k * db * db - restrwalls[ir].dMorse;
                dv = k * db / b;
            } else {
                const double fexp = exp(restrwalls[ir].aMorse * db);
                ener = restrwalls[ir].dMorse * (fexp * fexp - 2.0 * fexp);
                dv = -2.0 * restrwalls[ir].dMorse * restrwalls[ir].aMorse * (fexp - fexp * fexp) / b;
            }

            atomic_add_energy(&energy[E_RESTR_PRES], ener);

            atomic_add_force(&dvelocities[i].x, dv * dx);
            atomic_add_force(&dvelocities[i].y, dv * dy);
            atomic_add_force(&dvelocities[i].z, dv * dz);
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

    int blockSize = 256;
    int numBlocks = (host.n_restrwalls + blockSize - 1) / blockSize;
    calc_restrwall_forces_kernel<<<numBlocks, blockSize>>>(
        d_restrwalls,
        host.n_restrwalls,
        d_coords,
        host.energy.device(),
        d_dvelocities, d_heavy, host.topo);
}

void init_restrwall_force_kernel_data() {
}

void cleanup_restrwall_force() {
}
