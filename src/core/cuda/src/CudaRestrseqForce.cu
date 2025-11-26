#include "cuda/include/CudaRestrseqForce.cuh"
#include "utils.h"
#include "iostream"

namespace CudaRestrseqForce {
bool is_initialized = false;
restrseq_t* d_restrseq = nullptr;
coord_t* d_coords = nullptr;
coord_t* d_coords_top = nullptr;
double* d_upres_energy = nullptr;
atype_t* d_atypes = nullptr;
catype_t* d_catypes = nullptr;
bool* d_heavy = nullptr;
dvel_t* d_velocities = nullptr;

}  // namespace CudaRestrseqForce

__global__ void calc_restrseq_forces_kernel(
    int n_restrseqs,
    restrseq_t* restrseqs,
    coord_t* coords,
    coord_t* coords_top,
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
                dr.x += (coords[i].x - coords_top[i].x);
                dr.y += (coords[i].y - coords_top[i].y);
                dr.z += (coords[i].z - coords_top[i].z);
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
            if (heavy[i] || restrseqs[i].ih) {
                mass = catypes[atypes[i].code - 1].m;
                totmass += mass;
                dr.x += (coords[i].x - coords_top[i].x) * mass;
                dr.y += (coords[i].y - coords_top[i].y) * mass;
                dr.z += (coords[i].z - coords_top[i].z) * mass;
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
                dr.x = coords[i].x - coords_top[i].x;
                dr.y = coords[i].y - coords_top[i].y;
                dr.z = coords[i].z - coords_top[i].z;

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
    using namespace CudaRestrseqForce;
    if (!is_initialized) {
        check_cudaMalloc((void**)&d_restrseq, sizeof(restrseq_t) * n_restrseqs);
        check_cudaMalloc((void**)&d_coords, sizeof(coord_t) * n_atoms);
        check_cudaMalloc((void**)&d_coords_top, sizeof(coord_t) * n_atoms);
        check_cudaMalloc((void**)&d_upres_energy, sizeof(double));
        check_cudaMalloc((void**)&d_atypes, sizeof(atype_t) * n_atypes);
        check_cudaMalloc((void**)&d_catypes, sizeof(catype_t) * n_catypes);
        check_cudaMalloc((void**)&d_heavy, sizeof(bool) * n_atoms);
        check_cudaMalloc((void**)&d_velocities, sizeof(dvel_t) * n_atoms);

        cudaMemcpy(d_restrseq, restrseqs, sizeof(restrseq_t) * n_restrseqs, cudaMemcpyHostToDevice);
        cudaMemcpy(d_atypes, atypes, sizeof(atype_t) * n_atypes, cudaMemcpyHostToDevice);
        cudaMemcpy(d_catypes, catypes, sizeof(catype_t) * n_catypes, cudaMemcpyHostToDevice);
        cudaMemcpy(d_heavy, heavy, sizeof(bool) * n_atoms, cudaMemcpyHostToDevice);

        is_initialized = true;
    }

    cudaMemcpy(d_coords, coords, sizeof(coord_t) * n_atoms, cudaMemcpyHostToDevice);
    cudaMemcpy(d_coords_top, coords_top, sizeof(coord_t) * n_atoms, cudaMemcpyHostToDevice);
    cudaMemcpy(d_velocities, dvelocities, sizeof(dvel_t) * n_atoms, cudaMemcpyHostToDevice);
    cudaMemset(d_upres_energy, 0, sizeof(double));
    int blockSize = 256;
    int numBlocks = (n_restrseqs + blockSize - 1) / blockSize;
    calc_restrseq_forces_kernel<<<numBlocks, blockSize>>>(
        n_restrseqs,
        d_restrseq,
        d_coords,
        d_coords_top,
        d_atypes,
        d_catypes,
        d_heavy,
        d_velocities,
        d_upres_energy
    );
    cudaDeviceSynchronize();
    double upres_energy;
    cudaMemcpy(&upres_energy, d_upres_energy, sizeof(double), cudaMemcpyDeviceToHost);
    E_restraint.Upres = upres_energy;
    printf("Restrseq U_upres: %f\n", upres_energy);
    cudaMemcpy(dvelocities, d_velocities, sizeof(dvel_t) * n_atoms, cudaMemcpyDeviceToHost);
}

void cleanup_restrseq_force() {
    using namespace CudaRestrseqForce;
    if (is_initialized) {
        cudaFree(d_restrseq);
        cudaFree(d_coords);
        cudaFree(d_coords_top);
        cudaFree(d_upres_energy);
        cudaFree(d_atypes);
        cudaFree(d_catypes);
        cudaFree(d_heavy);
        cudaFree(d_velocities);
        is_initialized = false;
    }
}
