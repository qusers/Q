#include <iostream>
#include <math.h>
#include <vector>

#include "common/include/context.h"
#include "common/include/constants.h"
#include "cuda/include/cuda_force_accum.cuh"
#include "cuda/include/cuda_polx_water_force.cuh"
#include "cuda_utility.cuh"

namespace CudaPolxWaterForce {
bool is_initialized = false;

// in host
int* water_shell = nullptr;
int* water_rank = nullptr;
int* polx_list_sh = nullptr;  // use 1d array to simulate 2d array

real_t* d_energy;
int* d_list_sh = nullptr;
real_t* d_theta = nullptr;
real_t* d_theta0 = nullptr;
real_t* d_tdum = nullptr;
force_fixed_storage_t* d_avtheta_fixed = nullptr;
int* d_water_shell = nullptr;
int* d_water_rank = nullptr;

}  // namespace CudaPolxWaterForce

__global__ void calc_polx_theta_and_shells(
    int n_waters, int n_shells, int n_atoms_solute,
    coord_t* coords, topo_t topo, shell_t* wshells, int* list_sh,
    real_t* theta, real_t* theta0, real_t* tdum) {
    int idx = blockIdx.x * blockDim.x + threadIdx.x;
    if (idx >= n_waters) return;
    int i = idx;

    int wi, iis;
    coord_t rmu, rcu;
    real_t rm, rc;
    real_t cos_th;

    theta[i] = 0;
    theta0[i] = 0;

    wi = n_atoms_solute + 3 * i;

    rmu.x = coords[wi + 1].x + coords[wi + 2].x - 2 * coords[wi].x;
    rmu.y = coords[wi + 1].y + coords[wi + 2].y - 2 * coords[wi].y;
    rmu.z = coords[wi + 1].z + coords[wi + 2].z - 2 * coords[wi].z;

    rm = sqrt(rmu.x * rmu.x + rmu.y * rmu.y + rmu.z * rmu.z);

    rmu.x /= rm;
    rmu.y /= rm;
    rmu.z /= rm;

    rcu.x = coords[wi].x - topo.solvent_center.x;
    rcu.y = coords[wi].y - topo.solvent_center.y;
    rcu.z = coords[wi].z - topo.solvent_center.z;
    rc = sqrt(rcu.x * rcu.x + rcu.y * rcu.y + rcu.z * rcu.z);
    rcu.x /= rc;
    rcu.y /= rc;
    rcu.z /= rc;

    cos_th = rmu.x * rcu.x + rmu.y * rcu.y + rmu.z * rcu.z;
    if (cos_th > 1) cos_th = 1;
    if (cos_th < -1) cos_th = -1;
    theta[i] = acos(cos_th);
    tdum[i] = theta[i];
    // For waters outside inner shell, locate shell they're in
    if (rc > wshells[n_shells - 1].router - wshells[n_shells - 1].dr) {
        for (iis = n_shells - 1; iis > 0; iis--) {
            if (rc <= wshells[iis].router) break;
        }

        int pos = atomicAdd(&wshells[iis].n_inshell, 1);
        // printf("Water %d assigned to shell %d at position %d\n", i, iis, pos);
        // [pos][iis] -> (pos * n_shells + iis)
        list_sh[pos * n_shells + iis] = i;
    }
}

__global__ void calc_polx_water_forces_kernel(
    int n_waters, int n_atoms_solute, shell_t* wshells,
    coord_t* coords, cuda_dvel_t* dvelocities, topo_t topo,
    real_t* theta, md_t md, real_t* energy,
    int* water_rank, int* water_shell,
    force_fixed_storage_t* avtheta_fixed) {
    int idx = blockIdx.x * blockDim.x + threadIdx.x;
    if (idx >= n_waters) return;

    int il = water_rank[idx];
    int is = water_shell[idx];
    if (is < 0) return;  // water not in any shell this step

    int wi, ii;
    coord_t rmu, rcu, f1O, f1H1, f1H2, f2;
    real_t rm, rc;
    real_t cos_th;
    real_t avtdum, arg, f0, dv;
    real_t ener;

    avtdum = 0;
    ii = idx;
    arg = 1 + ((1 - 2 * (real_t)(il + 1)) / (real_t)wshells[is].n_inshell);
    real_t theta_val = acos(arg);
    theta_val = theta_val - 3 * sin(theta_val) * wshells[is].cstb / 2;
    if (theta_val < 0) theta_val = 0;
    if (theta_val > M_PI) theta_val = M_PI;

    avtdum += theta[ii];
    const real_t dtheta = theta[ii] - theta_val + wshells[is].theta_corr;
    ener = .5 * md.polarisation_force * dtheta * dtheta;
    // E_restraint.Upolx += ener;
    atomicAdd(energy, ener);

    dv = md.polarisation_force * dtheta;
    wi = n_atoms_solute + 3 * ii;

    rmu.x = coords[wi + 1].x + coords[wi + 2].x - 2 * coords[wi].x;
    rmu.y = coords[wi + 1].y + coords[wi + 2].y - 2 * coords[wi].y;
    rmu.z = coords[wi + 1].z + coords[wi + 2].z - 2 * coords[wi].z;

    rm = sqrt(rmu.x * rmu.x + rmu.y * rmu.y + rmu.z * rmu.z);

    rmu.x /= rm;
    rmu.y /= rm;
    rmu.z /= rm;

    rcu.x = coords[wi].x - topo.solvent_center.x;
    rcu.y = coords[wi].y - topo.solvent_center.y;
    rcu.z = coords[wi].z - topo.solvent_center.z;
    rc = sqrt(rcu.x * rcu.x + rcu.y * rcu.y + rcu.z * rcu.z);
    rcu.x /= rc;
    rcu.y /= rc;
    rcu.z /= rc;

    cos_th = rmu.x * rcu.x + rmu.y * rcu.y + rmu.z * rcu.z;
    if (cos_th > 1) cos_th = 1;
    if (cos_th < -1) cos_th = -1;
    f0 = sin(acos(cos_th));
    if (abs(f0) < k_singular_sin_epsilon) f0 = k_singular_sin_epsilon;
    f0 = -1.0 / f0;
    f0 *= dv;

    f1O.x = -2 * (rcu.x - rmu.x * cos_th) / rm;
    f1O.y = -2 * (rcu.y - rmu.y * cos_th) / rm;
    f1O.z = -2 * (rcu.z - rmu.z * cos_th) / rm;
    f1H1.x = (rcu.x - rmu.x * cos_th) / rm;
    f1H1.y = (rcu.y - rmu.y * cos_th) / rm;
    f1H1.z = (rcu.z - rmu.z * cos_th) / rm;
    f1H2.x = (rcu.x - rmu.x * cos_th) / rm;
    f1H2.y = (rcu.y - rmu.y * cos_th) / rm;
    f1H2.z = (rcu.z - rmu.z * cos_th) / rm;

    f2.x = (rmu.x - rcu.x * cos_th) / rc;
    f2.y = (rmu.y - rcu.y * cos_th) / rc;
    f2.z = (rmu.z - rcu.z * cos_th) / rc;

    atomic_add_force_component(&dvelocities[wi].x, f0 * (f1O.x + f2.x));
    atomic_add_force_component(&dvelocities[wi].y, f0 * (f1O.y + f2.y));
    atomic_add_force_component(&dvelocities[wi].z, f0 * (f1O.z + f2.z));
    atomic_add_force_component(&dvelocities[wi + 1].x, f0 * (f1H1.x));
    atomic_add_force_component(&dvelocities[wi + 1].y, f0 * (f1H1.y));
    atomic_add_force_component(&dvelocities[wi + 1].z, f0 * (f1H1.z));
    atomic_add_force_component(&dvelocities[wi + 2].x, f0 * (f1H2.x));
    atomic_add_force_component(&dvelocities[wi + 2].y, f0 * (f1H2.y));
    atomic_add_force_component(&dvelocities[wi + 2].z, f0 * (f1H2.z));

    atomic_add_force_component(&avtheta_fixed[is], avtdum / (real_t)wshells[is].n_inshell);
    atomicAdd(&wshells[is].avn_inshell, real_t(1));
}

void sort_waters() {
    using namespace CudaPolxWaterForce;
    auto& ctx = Context::instance();
    auto *wshells = ctx.wshells->cpu_data_p;

    int imin, jmin, jw;
    real_t tmin;
    // Sort the waters according to theta
    for (int is = 0; is < ctx.n_shells; is++) {
        imin = 0;
        for (int il = 0; il < wshells[is].n_inshell; il++) {
            tmin = 2 * M_PI;
            for (int jl = 0; jl < wshells[is].n_inshell; jl++) {
                // printf("Searching water %d in shell %d, total number: %d\n", jl, is, wshells[is].n_inshell);
                jw = polx_list_sh[jl * ctx.n_shells + is];
                // printf("Water %d in shell %d has theta %f\n", jw, is, ctx.tdum[jw] * 180 / M_PI);
                if (ctx.tdum[jw] < tmin) {
                    jmin = jw;
                    tmin = ctx.tdum[jw];
                }
            }
            ctx.nsort[imin][is] = jmin;
            water_rank[jmin] = imin;
            water_shell[jmin] = is;
            imin++;
            ctx.tdum[jmin] = 99999;
        }
    }
}

void calc_polx_water_forces_host(int iteration) {
    auto& ctx = Context::instance();
    auto *wshells = ctx.wshells->cpu_data_p;

    using namespace CudaPolxWaterForce;
    if (ctx.n_shells > 0) {
        std::vector<force_fixed_storage_t> avtheta_fixed(ctx.n_shells);
        cudaMemcpy(avtheta_fixed.data(), d_avtheta_fixed, ctx.n_shells * sizeof(force_fixed_storage_t), cudaMemcpyDeviceToHost);
        for (int is = 0; is < ctx.n_shells; is++) {
            wshells[is].avtheta = fixed_to_force(force_fixed_from_storage(avtheta_fixed[is]));
        }
    }

    for (int is = 0; is < ctx.n_shells; is++) {
        wshells[is].n_inshell = 0;
        if (iteration == 0 && !ctx.has_restart_wshell_theta_corr) {
            wshells[is].theta_corr = 0;
        }
    }

    coord_t* d_coords = ctx.coords->gpu_data_p;
    auto d_dvelocities = cuda_force_accum_buffer(ctx);
    shell_t* d_wshells = ctx.wshells->gpu_data_p;
    ctx.wshells->upload();

    int blockSize = 256;
    int numBlocks = (ctx.n_waters + blockSize - 1) / blockSize;
    // printf("Calculated theta for %d waters in %d shells\n", ctx.n_waters, ctx.n_shells);
    calc_polx_theta_and_shells<<<numBlocks, blockSize>>>(
        ctx.n_waters, ctx.n_shells, ctx.n_atoms_solute, d_coords, ctx.topo, d_wshells, d_list_sh, d_theta, d_theta0, d_tdum);
    // printf("Calculated theta for %d waters in %d shells\n", ctx.n_waters, ctx.n_shells);

    // todo: sort in cpu now..
    ctx.wshells->download();
    cudaMemcpy(ctx.tdum.data(), d_tdum, ctx.n_waters * sizeof(real_t), cudaMemcpyDeviceToHost);
    cudaMemcpy(polx_list_sh, d_list_sh, ctx.n_max_inshell * ctx.n_shells * sizeof(int), cudaMemcpyDeviceToHost);

    // Reset per-water metadata; only waters placed in shells will be overwritten in sort_waters().
    for (int i = 0; i < ctx.n_waters; i++) {
        water_rank[i] = -1;
        water_shell[i] = -1;
    }

    // printf("Sorting %d waters in %d shells\n", ctx.n_waters, ctx.n_shells);
    sort_waters();
    // printf("Sorted %d waters in %d shells\n", ctx.n_waters, ctx.n_shells);

    cudaMemcpy(d_water_rank, water_rank, ctx.n_waters * sizeof(int), cudaMemcpyHostToDevice);
    cudaMemcpy(d_water_shell, water_shell, ctx.n_waters * sizeof(int), cudaMemcpyHostToDevice);

    // Update theta_corr, averages
    if (iteration != 0 && iteration % itdis_update == 0) {
        for (int is = 0; is < ctx.n_shells; is++) {
            printf("SHELL %d\n", is);
            wshells[is].avtheta /= (real_t)itdis_update;
            wshells[is].avn_inshell /= (real_t)itdis_update;
            wshells[is].theta_corr = wshells[is].theta_corr + wshells[is].avtheta - acos(wshells[is].cstb);
            printf("average theta = %f, average in shell = %f, theta_corr = %f\n",
                   wshells[is].avtheta * 180 / M_PI, wshells[is].avn_inshell, wshells[is].theta_corr * 180 / M_PI);
            wshells[is].avtheta = 0;
            wshells[is].avn_inshell = 0;
        }
        cudaMemset(d_avtheta_fixed, 0, ctx.n_shells * sizeof(force_fixed_storage_t));
        ctx.wshells->upload();
    }

    // Calculate energy and force
    cudaMemset(d_energy, 0, sizeof(real_t));
    calc_polx_water_forces_kernel<<<numBlocks, blockSize>>>(
        ctx.n_waters, ctx.n_atoms_solute, d_wshells, d_coords, d_dvelocities, ctx.topo,
        d_theta, ctx.md, d_energy, d_water_rank, d_water_shell, d_avtheta_fixed);
    real_t energy;
    cudaMemcpy(&energy, d_energy, sizeof(real_t), cudaMemcpyDeviceToHost);
    ctx.E_restraint.Upolx += energy;
    ctx.wshells->download();
    if (ctx.n_shells > 0) {
        std::vector<force_fixed_storage_t> avtheta_fixed(ctx.n_shells);
        cudaMemcpy(avtheta_fixed.data(), d_avtheta_fixed, ctx.n_shells * sizeof(force_fixed_storage_t), cudaMemcpyDeviceToHost);
        for (int is = 0; is < ctx.n_shells; is++) {
            wshells[is].avtheta = fixed_to_force(force_fixed_from_storage(avtheta_fixed[is]));
        }
    }
    // Copy back forces for all atoms (solute + solvent); water forces were being dropped.
}

void init_polx_water_force_kernel_data() {
    using namespace CudaPolxWaterForce;
    auto& ctx = Context::instance();
    if (!is_initialized) {
        water_rank = new int[ctx.n_waters];
        water_shell = new int[ctx.n_waters];
        polx_list_sh = new int[ctx.n_max_inshell * ctx.n_shells];

        check_cudaMalloc((void**)&d_energy, sizeof(real_t));
        check_cudaMalloc((void**)&d_list_sh, ctx.n_max_inshell * ctx.n_shells * sizeof(int));
        check_cudaMalloc((void**)&d_theta, ctx.n_waters * sizeof(real_t));
        check_cudaMalloc((void**)&d_theta0, ctx.n_waters * sizeof(real_t));
        check_cudaMalloc((void**)&d_tdum, ctx.n_waters * sizeof(real_t));
        check_cudaMalloc((void**)&d_avtheta_fixed, ctx.n_shells * sizeof(force_fixed_storage_t));
        cudaMemset(d_avtheta_fixed, 0, ctx.n_shells * sizeof(force_fixed_storage_t));
        check_cudaMalloc((void**)&d_water_rank, ctx.n_waters * sizeof(int));
        check_cudaMalloc((void**)&d_water_shell, ctx.n_waters * sizeof(int));

        is_initialized = true;
    }
}

void cleanup_polx_water_force() {
    using namespace CudaPolxWaterForce;
    if (is_initialized) {
        delete[] water_rank;
        delete[] water_shell;
        delete[] polx_list_sh;

        cudaFree(d_energy);
        cudaFree(d_list_sh);
        cudaFree(d_theta);
        cudaFree(d_theta0);
        cudaFree(d_tdum);
        cudaFree(d_avtheta_fixed);
        cudaFree(d_water_rank);
        cudaFree(d_water_shell);
        is_initialized = false;
    }
}
