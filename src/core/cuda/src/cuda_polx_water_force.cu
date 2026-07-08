#include <iostream>
#include <math.h>

#include "common/include/context.h"
#include "common/include/constants.h"
#include "cuda/include/cuda_polx_water_force.cuh"
#include "cuda_force_accumulation.cuh"

namespace CudaPolxWaterForce {
bool is_initialized = false;

// in host
int* water_shell = nullptr;
int* water_rank = nullptr;
int* polx_list_sh = nullptr;  // use 1d array to simulate 2d array

int* d_list_sh = nullptr;
double* d_theta = nullptr;
double* d_theta0 = nullptr;
double* d_tdum = nullptr;
int* d_water_shell = nullptr;
int* d_water_rank = nullptr;

}  // namespace CudaPolxWaterForce

__global__ void calc_polx_theta_and_shells(
    int n_waters, int n_shells, int n_atoms_solute,
    coord_t* coords, topo_t topo, shell_t* wshells, int* list_sh,
    double* theta, double* theta0, double* tdum) {
    int idx = blockIdx.x * blockDim.x + threadIdx.x;
    if (idx >= n_waters) return;
    int i = idx;

    int wi, iis;

    theta[i] = 0.0;
    theta0[i] = 0.0;

    wi = n_atoms_solute + 3 * i;

    double rmu_x = coords[wi + 1].x +
                   coords[wi + 2].x -
                   2.0 * coords[wi].x;
    double rmu_y = coords[wi + 1].y +
                   coords[wi + 2].y -
                   2.0 * coords[wi].y;
    double rmu_z = coords[wi + 1].z +
                   coords[wi + 2].z -
                   2.0 * coords[wi].z;

    const double rm = sqrt(rmu_x * rmu_x + rmu_y * rmu_y + rmu_z * rmu_z);

    rmu_x /= rm;
    rmu_y /= rm;
    rmu_z /= rm;

    double rcu_x = coords[wi].x - topo.solvent_center.x;
    double rcu_y = coords[wi].y - topo.solvent_center.y;
    double rcu_z = coords[wi].z - topo.solvent_center.z;
    const double rc = sqrt(rcu_x * rcu_x + rcu_y * rcu_y + rcu_z * rcu_z);
    rcu_x /= rc;
    rcu_y /= rc;
    rcu_z /= rc;

    double cos_th = rmu_x * rcu_x + rmu_y * rcu_y + rmu_z * rcu_z;
    if (cos_th > 1.0) cos_th = 1.0;
    if (cos_th < -1.0) cos_th = -1.0;
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
    coord_t* coords, dvel_t* dvelocities, topo_t topo,
    double* theta, md_t md, energy_accum_t* energy,
    int* water_rank, int* water_shell) {
    int idx = blockIdx.x * blockDim.x + threadIdx.x;
    if (idx >= n_waters) return;

    int il = water_rank[idx];
    int is = water_shell[idx];
    if (is < 0) return;  // water not in any shell this step

    int wi, ii;

    double avtdum = 0.0;
    ii = idx;
    const double arg = 1.0 + ((1.0 - 2.0 * static_cast<double>(il + 1)) /
                              static_cast<double>(wshells[is].n_inshell));
    double theta_val = acos(arg);
    theta_val = theta_val - 3.0 * sin(theta_val) * wshells[is].cstb / 2.0;
    if (theta_val < 0.0) theta_val = 0.0;
    if (theta_val > M_PI) theta_val = M_PI;

    avtdum += theta[ii];
    const double dtheta = theta[ii] - theta_val + wshells[is].theta_corr;
    const double ener = 0.5 * md.polarisation_force * dtheta * dtheta;
    // E_restraint.Upolx += ener;
    atomic_add_energy(&energy[E_RESTR_POLX], ener);

    const double dv = md.polarisation_force * dtheta;
    wi = n_atoms_solute + 3 * ii;

    double rmu_x = coords[wi + 1].x +
                   coords[wi + 2].x -
                   2.0 * coords[wi].x;
    double rmu_y = coords[wi + 1].y +
                   coords[wi + 2].y -
                   2.0 * coords[wi].y;
    double rmu_z = coords[wi + 1].z +
                   coords[wi + 2].z -
                   2.0 * coords[wi].z;

    const double rm = sqrt(rmu_x * rmu_x + rmu_y * rmu_y + rmu_z * rmu_z);

    rmu_x /= rm;
    rmu_y /= rm;
    rmu_z /= rm;

    double rcu_x = coords[wi].x - topo.solvent_center.x;
    double rcu_y = coords[wi].y - topo.solvent_center.y;
    double rcu_z = coords[wi].z - topo.solvent_center.z;
    const double rc = sqrt(rcu_x * rcu_x + rcu_y * rcu_y + rcu_z * rcu_z);
    rcu_x /= rc;
    rcu_y /= rc;
    rcu_z /= rc;

    double cos_th = rmu_x * rcu_x + rmu_y * rcu_y + rmu_z * rcu_z;
    if (cos_th > 1.0) cos_th = 1.0;
    if (cos_th < -1.0) cos_th = -1.0;
    double f0 = sin(acos(cos_th));
    const double sin_epsilon = k_singular_sin_epsilon;
    if (fabs(f0) < sin_epsilon) f0 = sin_epsilon;
    f0 = -1.0 / f0;
    f0 *= dv;

    const double f1O_x = -2.0 * (rcu_x - rmu_x * cos_th) / rm;
    const double f1O_y = -2.0 * (rcu_y - rmu_y * cos_th) / rm;
    const double f1O_z = -2.0 * (rcu_z - rmu_z * cos_th) / rm;
    const double f1H1_x = (rcu_x - rmu_x * cos_th) / rm;
    const double f1H1_y = (rcu_y - rmu_y * cos_th) / rm;
    const double f1H1_z = (rcu_z - rmu_z * cos_th) / rm;
    const double f1H2_x = (rcu_x - rmu_x * cos_th) / rm;
    const double f1H2_y = (rcu_y - rmu_y * cos_th) / rm;
    const double f1H2_z = (rcu_z - rmu_z * cos_th) / rm;

    const double f2_x = (rmu_x - rcu_x * cos_th) / rc;
    const double f2_y = (rmu_y - rcu_y * cos_th) / rc;
    const double f2_z = (rmu_z - rcu_z * cos_th) / rc;

    atomic_add_force(&dvelocities[wi].x, f0 * (f1O_x + f2_x));
    atomic_add_force(&dvelocities[wi].y, f0 * (f1O_y + f2_y));
    atomic_add_force(&dvelocities[wi].z, f0 * (f1O_z + f2_z));
    atomic_add_force(&dvelocities[wi + 1].x, f0 * f1H1_x);
    atomic_add_force(&dvelocities[wi + 1].y, f0 * f1H1_y);
    atomic_add_force(&dvelocities[wi + 1].z, f0 * f1H1_z);
    atomic_add_force(&dvelocities[wi + 2].x, f0 * f1H2_x);
    atomic_add_force(&dvelocities[wi + 2].y, f0 * f1H2_y);
    atomic_add_force(&dvelocities[wi + 2].z, f0 * f1H2_z);

    atomicAdd(&wshells[is].avtheta, avtdum / static_cast<double>(wshells[is].n_inshell));
    if (il == 0) {
        atomicAdd(&wshells[is].avn_inshell, static_cast<double>(wshells[is].n_inshell));
    }
}

void sort_waters() {
    using namespace CudaPolxWaterForce;
    auto& ctx = Context::instance();
    auto *wshells = ctx.wshells->cpu_data_p;

    int imin, jmin, jw;
    double tmin;
    // Sort the waters according to theta
    for (int is = 0; is < ctx.n_shells; is++) {
        imin = 0;
        for (int il = 0; il < wshells[is].n_inshell; il++) {
            tmin = 2.0 * M_PI;
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
            ctx.tdum[jmin] = 99999.0;
        }
    }
}

void calc_polx_water_forces_host(int iteration) {
    auto& ctx = Context::instance();
    auto *wshells = ctx.wshells->cpu_data_p;

    for (int is = 0; is < ctx.n_shells; is++) {
        wshells[is].n_inshell = 0;
    }
    using namespace CudaPolxWaterForce;

    coord_t* d_coords = ctx.coords->gpu_data_p;
    dvel_t* d_dvelocities = ctx.dvelocities->gpu_data_p;
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
    cudaMemcpy(ctx.tdum.data(), d_tdum, ctx.n_waters * sizeof(double), cudaMemcpyDeviceToHost);
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
            wshells[is].avtheta /= static_cast<double>(itdis_update);
            wshells[is].avn_inshell /= static_cast<double>(itdis_update);
            wshells[is].theta_corr = wshells[is].theta_corr + wshells[is].avtheta - acos(wshells[is].cstb);
            printf("average theta = %f, average in shell = %f, theta_corr = %f\n",
                   wshells[is].avtheta * 180 / M_PI, wshells[is].avn_inshell, wshells[is].theta_corr * 180 / M_PI);
            wshells[is].avtheta = 0.0;
            wshells[is].avn_inshell = 0.0;
        }
        ctx.wshells->upload();
    }

    // Calculate energy and force
    calc_polx_water_forces_kernel<<<numBlocks, blockSize>>>(
        ctx.n_waters, ctx.n_atoms_solute, d_wshells, d_coords, d_dvelocities, ctx.topo,
        d_theta, ctx.md, ctx.energy.device(), d_water_rank, d_water_shell);
    ctx.wshells->download();
    // Copy back forces for all atoms (solute + solvent); water forces were being dropped.
}

void init_polx_water_force_kernel_data() {
    using namespace CudaPolxWaterForce;
    auto& ctx = Context::instance();
    if (!is_initialized) {
        water_rank = new int[ctx.n_waters];
        water_shell = new int[ctx.n_waters];
        polx_list_sh = new int[ctx.n_max_inshell * ctx.n_shells];

        check_cudaMalloc((void**)&d_list_sh, ctx.n_max_inshell * ctx.n_shells * sizeof(int));
        check_cudaMalloc((void**)&d_theta, ctx.n_waters * sizeof(double));
        check_cudaMalloc((void**)&d_theta0, ctx.n_waters * sizeof(double));
        check_cudaMalloc((void**)&d_tdum, ctx.n_waters * sizeof(double));
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

        cudaFree(d_list_sh);
        cudaFree(d_theta);
        cudaFree(d_theta0);
        cudaFree(d_tdum);
        cudaFree(d_water_rank);
        cudaFree(d_water_shell);
        is_initialized = false;
    }
}
