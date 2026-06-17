#include "cuda/include/cuda_restrseq_force.cuh"
#include "common/include/context.h"
#include "iostream"
#include "cuda_force_accumulation.cuh"


namespace CudaRestrseqForce {
bool is_initialized = false;
energy_accum_t* d_upres_energy;
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
    energy_accum_t* upres_energy) {
    int s = blockIdx.x * blockDim.x + threadIdx.x;
    if (s >= n_restrseqs) return;

    const double k = restrseqs[s].k;

    double dx = 0.0;
    double dy = 0.0;
    double dz = 0.0;
    int n_ctr = 0;
    double totmass = 0.0;

    // Geometric center
    if (restrseqs[s].to_center == 1) {
        for (int i = restrseqs[s].ai - 1; i < restrseqs[s].aj - 1; i++) {
            if (heavy[i] || restrseqs[s].ih) {
                n_ctr++;
                dx += coords[i].x - coords_init[i].x;
                dy += coords[i].y - coords_init[i].y;
                dz += coords[i].z - coords_init[i].z;
            }
        }

        if (n_ctr > 0) {
            const double inv_n_ctr = 1.0 / static_cast<double>(n_ctr);
            dx *= inv_n_ctr;
            dy *= inv_n_ctr;
            dz *= inv_n_ctr;
            const double r2 = dx * dx + dy * dy + dz * dz;
            const double ener = 0.5 * k * r2;
            atomic_add_energy(upres_energy, ener);

            for (int i = restrseqs[s].ai - 1; i < restrseqs[s].aj - 1; i++) {
                if (heavy[i] || restrseqs[s].ih) {
                    const double mass = catypes[atypes[i].code - 1].m;
                    const double tmp = 12.010;
                    atomic_add_force(&dvelocities[i].x, k * dx * mass / tmp);
                    atomic_add_force(&dvelocities[i].y, k * dy * mass / tmp);
                    atomic_add_force(&dvelocities[i].z, k * dz * mass / tmp);
                }
            }
        }
    }

    // Mass center
    else if (restrseqs[s].to_center == 2) {
        for (int i = restrseqs[s].ai - 1; i < restrseqs[s].aj - 1; i++) {
            if (heavy[i] || restrseqs[s].ih) {
                const double mass = catypes[atypes[i].code - 1].m;
                totmass += mass;
                dx += (coords[i].x - coords_init[i].x) * mass;
                dy += (coords[i].y - coords_init[i].y) * mass;
                dz += (coords[i].z - coords_init[i].z) * mass;
            }
        }

        if (totmass > 0) {
            dx /= totmass;
            dy /= totmass;
            dz /= totmass;
            const double r2 = dx * dx + dy * dy + dz * dz;
            const double ener = 0.5 * k * r2;
            atomic_add_energy(upres_energy, ener);

            for (int i = restrseqs[s].ai - 1; i < restrseqs[s].aj - 1; i++) {
                if (heavy[i] || restrseqs[s].ih) {
                    atomic_add_force(&dvelocities[i].x, k * dx);
                    atomic_add_force(&dvelocities[i].y, k * dy);
                    atomic_add_force(&dvelocities[i].z, k * dz);
                }
            }
        }
    }

    // Restrain to topology coordinate
    else {
        for (int i = restrseqs[s].ai - 1; i < restrseqs[s].aj - 1; i++) {
            if (heavy[i] || restrseqs[s].ih) {
                const double dx_atom = coords[i].x - coords_init[i].x;
                const double dy_atom = coords[i].y - coords_init[i].y;
                const double dz_atom = coords[i].z - coords_init[i].z;

                const double r2 = dx_atom * dx_atom + dy_atom * dy_atom + dz_atom * dz_atom;
                const double ener = 0.5 * k * r2;
                atomic_add_energy(upres_energy, ener);

                atomic_add_force(&dvelocities[i].x, k * dx_atom);
                atomic_add_force(&dvelocities[i].y, k * dy_atom);
                atomic_add_force(&dvelocities[i].z, k * dz_atom);
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
    cudaMemset(d_upres_energy, 0, sizeof(energy_accum_t));
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
    energy_accum_t upres_energy;
    cudaMemcpy(&upres_energy, d_upres_energy, sizeof(energy_accum_t), cudaMemcpyDeviceToHost);
    const double energy = energy_from_accum(upres_energy);
    host.E_restraint.Upres = energy;
    printf("Restrseq U_upres: %f\n", energy);
}

void init_restrseq_force_kernel_data() {
    using namespace CudaRestrseqForce;
    if (!is_initialized) {
        check_cudaMalloc((void**)&d_upres_energy, sizeof(energy_accum_t));
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
