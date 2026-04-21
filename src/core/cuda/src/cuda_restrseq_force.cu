#include "cuda/include/cuda_restrseq_force.cuh"
#include "common/include/context.h"
#include "iostream"

namespace CudaRestrseqForce {
bool is_initialized = false;
double* d_upres_energy;
}  // namespace CudaRestrseqForce
__global__ void calc_restrseq_forces_kernel(
    int n_restrseqs,
    restrseq_t* restrseqs,
    coord_t* coords,
    coord_t* coords_init,
    atype_t* atypes,
    catype_t* catypes,
    bool* heavy,
    dvel_t* dvelocities,
    double* upres_energy) {
    int s = blockIdx.x * blockDim.x + threadIdx.x;
    if (s >= n_restrseqs) return;

    double k, mass, totmass;
    coord_t dr;
    double r2, ener;

    k = restrseqs[s].k;

    dr.x = 0;
    dr.y = 0;
    dr.z = 0;
    int n_ctr = 0;
    totmass = 0;

    // Geometric center
    if (restrseqs[s].to_center == 1) {
        for (int i = restrseqs[s].ai - 1; i < restrseqs[s].aj - 1; i++) {
            if (heavy[i] || restrseqs[s].ih) {
                n_ctr++;
                dr.x += (coords[i].x - coords_init[i].x);
                dr.y += (coords[i].y - coords_init[i].y);
                dr.z += (coords[i].z - coords_init[i].z);
            }
        }

        if (n_ctr > 0) {
            dr.x /= n_ctr;
            dr.y /= n_ctr;
            dr.z /= n_ctr;
            r2 = pow(dr.x, 2) + pow(dr.y, 2) + pow(dr.z, 2);
            ener = .5 * k * r2;
            atomicAdd(upres_energy, ener);

            for (int i = restrseqs[s].ai - 1; i < restrseqs[s].aj - 1; i++) {
                if (heavy[i] || restrseqs[s].ih) {
                    mass = catypes[atypes[i].code - 1].m;
                    atomicAdd(&dvelocities[i].x, (k * dr.x * mass / 12.010));
                    atomicAdd(&dvelocities[i].y, (k * dr.y * mass / 12.010));
                    atomicAdd(&dvelocities[i].z, (k * dr.z * mass / 12.010));
                }
            }
        }
    }

    // Mass center
    else if (restrseqs[s].to_center == 2) {
        for (int i = restrseqs[s].ai - 1; i < restrseqs[s].aj - 1; i++) {
            if (heavy[i] || restrseqs[s].ih) {
                mass = catypes[atypes[i].code - 1].m;
                totmass += mass;
                dr.x += (coords[i].x - coords_init[i].x) * mass;
                dr.y += (coords[i].y - coords_init[i].y) * mass;
                dr.z += (coords[i].z - coords_init[i].z) * mass;
            }
        }

        if (totmass > 0) {
            dr.x /= totmass;
            dr.y /= totmass;
            dr.z /= totmass;
            r2 = pow(dr.x, 2) + pow(dr.y, 2) + pow(dr.z, 2);
            ener = .5 * k * r2;
            atomicAdd(upres_energy, ener);

            for (int i = restrseqs[s].ai - 1; i < restrseqs[s].aj - 1; i++) {
                if (heavy[i] || restrseqs[s].ih) {
                    mass = catypes[atypes[i].code - 1].m;
                    atomicAdd(&dvelocities[i].x, k * dr.x);
                    atomicAdd(&dvelocities[i].y, k * dr.y);
                    atomicAdd(&dvelocities[i].z, k * dr.z);
                }
            }
        }
    }

    // Restrain to topology coordinate
    else {
        for (int i = restrseqs[s].ai - 1; i < restrseqs[s].aj - 1; i++) {
            if (heavy[i] || restrseqs[s].ih) {
                dr.x = coords[i].x - coords_init[i].x;
                dr.y = coords[i].y - coords_init[i].y;
                dr.z = coords[i].z - coords_init[i].z;

                r2 = pow(dr.x, 2) + pow(dr.y, 2) + pow(dr.z, 2);
                ener = .5 * k * r2;
                atomicAdd(upres_energy, ener);

                atomicAdd(&dvelocities[i].x, k * dr.x);
                atomicAdd(&dvelocities[i].y, k * dr.y);
                atomicAdd(&dvelocities[i].z, k * dr.z);
            }
        }
    }
}

void calc_restrseq_forces_host() {
    auto& host = Context::instance();
    if (host.n_restrseqs == 0) return;
    using namespace CudaRestrseqForce;
    auto d_restrseq = host.restrseqs->gpu_data_p;
    auto d_coords = host.coords->gpu_data_p;
    auto d_coords_init = host.coords_init->gpu_data_p;
    auto d_atypes = host.atypes->gpu_data_p;
    auto d_catypes = host.catypes->gpu_data_p;
    auto d_heavy = host.heavy->gpu_data_p;
    auto d_dvelocities = host.dvelocities->gpu_data_p;
    cudaMemset(d_upres_energy, 0, sizeof(double));
    // ctx.sync_all_to_device();

    int blockSize = 256;
    int numBlocks = (host.n_restrseqs + blockSize - 1) / blockSize;
    calc_restrseq_forces_kernel<<<numBlocks, blockSize>>>(
        host.n_restrseqs,
        d_restrseq,
        d_coords,
        d_coords_init,
        d_atypes,
        d_catypes,
        d_heavy,
        d_dvelocities,
        d_upres_energy);
    cudaDeviceSynchronize();
    double upres_energy;
    cudaMemcpy(&upres_energy, d_upres_energy, sizeof(double), cudaMemcpyDeviceToHost);
    host.E_restraint.Upres = upres_energy;
    printf("Restrseq U_upres: %f\n", upres_energy);
}

void init_restrseq_force_kernel_data() {
    using namespace CudaRestrseqForce;
    if (!is_initialized) {
        check_cudaMalloc((void**)&d_upres_energy, sizeof(double));
        is_initialized = true;
    }
}

void cleanup_restrseq_force() {
    using namespace CudaRestrseqForce;
    if (is_initialized) {
        cudaFree(d_upres_energy);
        is_initialized = false;
    }
}
