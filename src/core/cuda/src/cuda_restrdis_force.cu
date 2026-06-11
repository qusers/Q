#include <iostream>
#include <vector>

#include "cuda/include/cuda_restrdis_force.cuh"
#include "cuda/include/cuda_utility.cuh"
#include "common/include/context.h"
#include "cuda_force_accumulation.cuh"
namespace CudaRestrdisForce {
bool is_initialized = false;
energy_accum_t* d_urestr;
energy_accum_t* d_upres;
}  // namespace CudaRestrdisForce

__global__ void calc_restrdis_forces_kernel(
    restrdis_t* restrdists,
    int n_restrdists,
    coord_t* coords,
    real_t* lambdas,
    int n_lambdas,
    dvel_t* dvelocities,
    energy_accum_t* urestr,
    energy_accum_t* upres) {
    int idx = blockIdx.x * blockDim.x + threadIdx.x;
    if (idx >= n_restrdists) return;

    int state, i, j;

    int ir = idx;

    state = restrdists[ir].ipsi - 1;
    i = restrdists[ir].ai - 1;
    j = restrdists[ir].aj - 1;

    const double dx = static_cast<double>(coords[j].x) - static_cast<double>(coords[i].x);
    const double dy = static_cast<double>(coords[j].y) - static_cast<double>(coords[i].y);
    const double dz = static_cast<double>(coords[j].z) - static_cast<double>(coords[i].z);

    double lambda;
    if (restrdists[ir].ipsi != 0) {
        lambda = static_cast<double>(lambdas[state]);
    } else {
        lambda = 1.0;
    }

    const double b = sqrt(dx * dx + dy * dy + dz * dz);
    double db;
    if (b < restrdists[ir].d1) {
        db = b - restrdists[ir].d1;
    } else if (b > restrdists[ir].d2) {
        db = b - restrdists[ir].d2;
    } else {
        return;
    }

    const double ener = 0.5 * restrdists[ir].k * db * db;
    const double dv = lambda * restrdists[ir].k * db / b;

    atomic_add_force(&dvelocities[j].x, dx * dv);
    atomic_add_force(&dvelocities[j].y, dy * dv);
    atomic_add_force(&dvelocities[j].z, dz * dv);
    atomic_add_force(&dvelocities[i].x, -dx * dv);
    atomic_add_force(&dvelocities[i].y, -dy * dv);
    atomic_add_force(&dvelocities[i].z, -dz * dv);

    if (restrdists[ir].ipsi == 0) {
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

void calc_restrdis_forces_host() {
    auto& host = Context::instance();
    if (host.n_restrdists == 0) return;
    using namespace CudaRestrdisForce;
    auto d_restrdists = host.restrdists->gpu_data_p;
    auto d_coords = host.coords->gpu_data_p;
    auto d_lambdas = host.lambdas->gpu_data_p;
    auto d_dvelocities = host.dvelocities->gpu_data_p;

    if (host.n_lambdas > 0) {
        cudaMemset(d_urestr, 0, sizeof(energy_accum_t) * host.n_lambdas);
    }
    cudaMemset(d_upres, 0, sizeof(energy_accum_t));

    int blockSize = 256;
    int numBlocks = (host.n_restrdists + blockSize - 1) / blockSize;
    calc_restrdis_forces_kernel<<<numBlocks, blockSize>>>(
        d_restrdists,
        host.n_restrdists,
        d_coords,
        d_lambdas,
        host.n_lambdas,
        d_dvelocities,
        d_urestr,
        d_upres);
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

void init_restrdis_force_kernel_data() {
    using namespace CudaRestrdisForce;
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

void cleanup_restrdis_force() {
    using namespace CudaRestrdisForce;
    if (is_initialized) {
        if (d_urestr != nullptr) cudaFree(d_urestr);
        cudaFree(d_upres);
        d_urestr = nullptr;
        d_upres = nullptr;
        is_initialized = false;
    }
}
