#include <iostream>
#include <vector>

#include "cuda/include/cuda_restrpos_force.cuh"
#include "cuda/include/cuda_utility.cuh"
#include "common/include/context.h"
#include "cuda_force_accumulation.cuh"

namespace CudaRestrposForce {
bool is_initialized = false;
energy_accum_t* d_urestr;
energy_accum_t* d_upres;
}  // namespace CudaRestrposForce

__global__ void calc_restrpos_forces_kernel(
    restrpos_t* restrspos,
    int n_restrspos,
    coord_t* coords,
    real_t* lambdas,
    int n_lambdas,
    energy_accum_t* urestr,
    energy_accum_t* upres,
    dvel_t* dvelocities) {
    int idx = blockIdx.x * blockDim.x + threadIdx.x;
    if (idx >= n_restrspos) return;
    int ir = idx;

    int state, i;

    state = restrspos[ir].ipsi - 1;
    i = restrspos[ir].a - 1;

    const double dx = static_cast<double>(coords[i].x) - static_cast<double>(restrspos[ir].x.x);
    const double dy = static_cast<double>(coords[i].y) - static_cast<double>(restrspos[ir].x.y);
    const double dz = static_cast<double>(coords[i].z) - static_cast<double>(restrspos[ir].x.z);

    double lambda;
    if (restrspos[ir].ipsi != 0) {
        lambda = static_cast<double>(lambdas[state]);
    } else {
        lambda = 1.0;
    }

    const double kx = static_cast<double>(restrspos[ir].k.x);
    const double ky = static_cast<double>(restrspos[ir].k.y);
    const double kz = static_cast<double>(restrspos[ir].k.z);

    const double ener = 0.5 * kx * dx * dx + 0.5 * ky * dy * dy + 0.5 * kz * dz * dz;

    atomic_add_force(&dvelocities[i].x, kx * dx * lambda);
    atomic_add_force(&dvelocities[i].y, ky * dy * lambda);
    atomic_add_force(&dvelocities[i].z, kz * dz * lambda);

    if (restrspos[ir].ipsi == 0) {
        for (int k = 0; k < n_lambdas; k++) {
            atomic_add_energy(&urestr[k], ener);
        }
        if (n_lambdas == 0) {
            atomic_add_energy(upres, ener);
        }
    } else {
        atomic_add_energy(&urestr[state], ener);
    }
}
void calc_restrpos_forces_host() {
    auto& host = Context::instance();
    if (host.n_restrspos == 0) return;
    using namespace CudaRestrposForce;
    if (host.n_lambdas > 0) {
        cudaMemset(d_urestr, 0, sizeof(energy_accum_t) * host.n_lambdas);
    }
    cudaMemset(d_upres, 0, sizeof(energy_accum_t));

    auto d_restrspos = host.restrspos->gpu_data_p;
    auto d_coords = host.coords->gpu_data_p;
    auto d_lambdas = host.lambdas->gpu_data_p;
    auto d_dvelocities = host.dvelocities->gpu_data_p;

    int blockSize = 256;
    int numBlocks = (host.n_restrspos + blockSize - 1) / blockSize;
    calc_restrpos_forces_kernel<<<numBlocks, blockSize>>>(
        d_restrspos,
        host.n_restrspos,
        d_coords,
        d_lambdas,
        host.n_lambdas,
        d_urestr,
        d_upres,
        d_dvelocities);
    cudaDeviceSynchronize();
    if (host.n_lambdas > 0) {
        std::vector<energy_accum_t> urestr(host.n_lambdas, 0);
        cudaMemcpy(urestr.data(), d_urestr, sizeof(energy_accum_t) * host.n_lambdas, cudaMemcpyDeviceToHost);
        auto* EQ_restraint = host.EQ_restraint->cpu_data_p;
        for (int state = 0; state < host.n_lambdas; state++) {
            EQ_restraint[state].Urestr += energy_from_accum(urestr[state]);
        }
    }
    energy_accum_t upres = 0;
    cudaMemcpy(&upres, d_upres, sizeof(energy_accum_t), cudaMemcpyDeviceToHost);
    host.E_restraint.Upres += energy_from_accum(upres);
}

void init_restrpos_force_kernel_data() {
    using namespace CudaRestrposForce;
    if (!is_initialized) {
        auto& host = Context::instance();
        if (host.n_lambdas > 0) {
            check_cudaMalloc((void**)&d_urestr, sizeof(energy_accum_t) * host.n_lambdas);
        } else {
            d_urestr = nullptr;
        }
        check_cudaMalloc((void**)&d_upres, sizeof(energy_accum_t));
        is_initialized = true;
    }
}

void cleanup_restrpos_force() {
    using namespace CudaRestrposForce;
    if (is_initialized) {
        if (d_urestr != nullptr) cudaFree(d_urestr);
        cudaFree(d_upres);
        d_urestr = nullptr;
        d_upres = nullptr;
        is_initialized = false;
    }
}
