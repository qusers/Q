#include <vector>

#include "cuda/include/cuda_restrang_force.cuh"
#include "cuda/include/cuda_utility.cuh"
#include "common/include/context.h"
#include "cuda_force_accumulation.cuh"
namespace CudaRestrangForce {
bool is_initialized = false;
energy_accum_t* d_urestr;
energy_accum_t* d_upres;
}  // namespace CudaRestrangForce

__global__ void calc_restrang_force_kernel(
    restrang_t* restrangs,
    int n_restrangs,
    coord_t* coords,
    double* lambdas,
    int n_lambdas,
    dvel_t* dvelocities,
    energy_accum_t* urestr,
    energy_accum_t* upres) {
    int idx = blockIdx.x * blockDim.x + threadIdx.x;
    if (idx >= n_restrangs) return;
    int ir = idx;

    int state, i, j, k;

    state = restrangs[ir].ipsi - 1;
    i = restrangs[ir].ai - 1;
    j = restrangs[ir].aj - 1;
    k = restrangs[ir].ak - 1;

    const double dr_x = coords[i].x - coords[j].x;
    const double dr_y = coords[i].y - coords[j].y;
    const double dr_z = coords[i].z - coords[j].z;

    const double dr2_x = coords[k].x - coords[j].x;
    const double dr2_y = coords[k].y - coords[j].y;
    const double dr2_z = coords[k].z - coords[j].z;

    double lambda;
    if (restrangs[ir].ipsi != 0) {
        lambda = lambdas[state];
    } else {
        lambda = 1.0;
    }

    const double r2ij = dr_x * dr_x + dr_y * dr_y + dr_z * dr_z;
    const double r2jk = dr2_x * dr2_x + dr2_y * dr2_y + dr2_z * dr2_z;

    const double rij = sqrt(r2ij);
    const double rjk = sqrt(r2jk);

    double cos_th = dr_x * dr2_x + dr_y * dr2_y + dr_z * dr2_z;
    cos_th /= rij * rjk;

    if (cos_th > 1.0) cos_th = 1.0;
    if (cos_th < -1.0) cos_th = -1.0;

    const double th = acos(cos_th);
    const double dth = th - to_radians_device(restrangs[ir].ang);

    const double ener = 0.5 * restrangs[ir].k * dth * dth;
    const double dv = lambda * restrangs[ir].k * dth;

    double f1 = sin(th);
    const double sin_epsilon = k_singular_sin_epsilon;
    if (fabs(f1) < sin_epsilon) {
        f1 = -1.0 / sin_epsilon;
    } else {
        f1 = -1 / f1;
    }

    const double di_x = f1 * (dr2_x / (rij * rjk) - cos_th * dr_x / r2ij);
    const double di_y = f1 * (dr2_y / (rij * rjk) - cos_th * dr_y / r2ij);
    const double di_z = f1 * (dr2_z / (rij * rjk) - cos_th * dr_z / r2ij);
    const double dk_x = f1 * (dr_x / (rij * rjk) - cos_th * dr2_x / r2jk);
    const double dk_y = f1 * (dr_y / (rij * rjk) - cos_th * dr2_y / r2jk);
    const double dk_z = f1 * (dr_z / (rij * rjk) - cos_th * dr2_z / r2jk);

    atomic_add_force(&dvelocities[i].x, dv * di_x);
    atomic_add_force(&dvelocities[i].y, dv * di_y);
    atomic_add_force(&dvelocities[i].z, dv * di_z);
    atomic_add_force(&dvelocities[k].x, dv * dk_x);
    atomic_add_force(&dvelocities[k].y, dv * dk_y);
    atomic_add_force(&dvelocities[k].z, dv * dk_z);
    atomic_add_force(&dvelocities[j].x, -dv * (di_x + dk_x));
    atomic_add_force(&dvelocities[j].y, -dv * (di_y + dk_y));
    atomic_add_force(&dvelocities[j].z, -dv * (di_z + dk_z));

    if (restrangs[ir].ipsi == 0) {
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

void calc_restrang_force_host() {
    auto& host = Context::instance();
    if (host.n_restrangs == 0) return;
    using namespace CudaRestrangForce;

    auto d_restrangs = host.restrangs->gpu_data_p;
    auto d_coords = host.coords->gpu_data_p;
    auto d_lambdas = host.lambdas->gpu_data_p;
    auto d_dvelocities = host.dvelocities->gpu_data_p;

    if (host.n_lambdas > 0) {
        cudaMemset(d_urestr, 0, sizeof(energy_accum_t) * host.n_lambdas);
    }
    cudaMemset(d_upres, 0, sizeof(energy_accum_t));

    int blockSize = 256;
    int numBlocks = (host.n_restrangs + blockSize - 1) / blockSize;
    calc_restrang_force_kernel<<<numBlocks, blockSize>>>(
        d_restrangs,
        host.n_restrangs,
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

void init_restrang_force_kernel_data() {
    using namespace CudaRestrangForce;
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

void cleanup_restrang_force() {
    using namespace CudaRestrangForce;
    if (is_initialized) {
        if (d_urestr != nullptr) cudaFree(d_urestr);
        cudaFree(d_upres);
        d_urestr = nullptr;
        d_upres = nullptr;
        is_initialized = false;
    }
}
