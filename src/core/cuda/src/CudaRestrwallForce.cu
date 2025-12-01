#include <iostream>

#include "cuda/include/CudaContext.cuh"
#include "cuda/include/CudaRestrwallForce.cuh"
#include "utils.h"

__global__ void calc_restrwall_forces_kernel(
    restrwall_t* restrwalls,
    int n_restrwalls,
    coord_t* coords,
    double* energies,
    dvel_t* dvelocities,
    bool* heavy, topo_t topo) {
    int idx = blockIdx.x * blockDim.x + threadIdx.x;
    if (idx >= n_restrwalls) return;

    double k, b, db, ener, dv, fexp;
    coord_t dr;

    int ir = idx;
    for (int i = restrwalls[ir].ai - 1; i < restrwalls[ir].aj - 1; i++) {
        if (heavy[i] || restrwalls[ir].ih) {
            dr.x = coords[i].x - topo.solvent_center.x;
            dr.y = coords[i].y - topo.solvent_center.y;
            dr.x = coords[i].x - topo.solvent_center.x;

            b = sqrt(pow(dr.x, 2) + pow(dr.y, 2) + pow(dr.z, 2));
            db = b - restrwalls[ir].d;

            if (db > 0) {
                ener = .5 * k * pow(db, 2) - restrwalls[ir].dMorse;
                dv = k * db / b;
            } else {
                fexp = exp(restrwalls[ir].aMorse * db);
                ener = restrwalls[ir].dMorse * (fexp * fexp - 2 * fexp);
                dv = -2 * restrwalls[ir].dMorse * restrwalls[ir].aMorse * (fexp - fexp * fexp) / b;
            }

            atomicAdd(energies, ener);

            atomicAdd(&dvelocities[i].x, dv * dr.x);
            atomicAdd(&dvelocities[i].y, dv * dr.y);
            atomicAdd(&dvelocities[i].z, dv * dr.z);
        }
    }
}

void calc_restrwall_forces_host() {
    if (n_restrwalls == 0) return;
    CudaContext& ctx = CudaContext::instance();
    auto d_restrwalls = ctx.d_restrwalls;
    auto d_coords = ctx.d_coords;
    auto d_dvelocities = ctx.d_dvelocities;
    auto d_heavy = ctx.d_heavy;
    double* d_energies;
    check_cudaMalloc((void**)&d_energies, sizeof(double));
    cudaMemset(d_energies, 0, sizeof(double));

    int blockSize = 256;
    int numBlocks = (n_restrwalls + blockSize - 1) / blockSize;
    calc_restrwall_forces_kernel<<<numBlocks, blockSize>>>(
        d_restrwalls,
        n_restrwalls,
        d_coords,
        d_energies,
        d_dvelocities, d_heavy, topo);
    cudaDeviceSynchronize();
    double h_energy;
    cudaMemcpy(&h_energy, d_energies, sizeof(double), cudaMemcpyDeviceToHost);
    cudaMemcpy(dvelocities, d_dvelocities, sizeof(dvel_t) * n_atoms, cudaMemcpyDeviceToHost);
    printf("Restrwall energy: %f\n", h_energy);
    E_restraint.Upres += h_energy;
    cudaFree(d_energies);
}