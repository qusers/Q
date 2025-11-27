#include <iostream>

#include "cuda/include/CudaContext.cuh"
#include "cuda/include/CudaPolxWaterForce.cuh"
#include "utils.h"

namespace CudaPolxWaterForce {
bool is_initialized = false;

// in host
int* water_shell = nullptr;
int* water_rank = nullptr;
int* polx_list_sh = nullptr;  // use 1d array to simulate 2d array

}  // namespace CudaPolxWaterForce

__global__ void calc_polx_theta_and_shells(
    int n_waters, int n_shells, int n_atoms_solute,
    coord_t* coords, topo_t topo, shell_t* wshells, int* list_sh,
    double* theta, double* theta0, double* tdum) {
    int idx = blockIdx.x * blockDim.x + threadIdx.x;
    if (idx >= n_waters) return;
    int i = idx;

    int wi, iis;
    coord_t rmu, rcu;
    double rm, rc;
    double cos_th;

    theta[i] = 0;
    theta0[i] = 0;

    wi = n_atoms_solute + 3 * i;

    rmu.x = coords[wi + 1].x + coords[wi + 2].x - 2 * coords[wi].x;
    rmu.y = coords[wi + 1].y + coords[wi + 2].y - 2 * coords[wi].y;
    rmu.z = coords[wi + 1].z + coords[wi + 2].z - 2 * coords[wi].z;

    rm = sqrt(pow(rmu.x, 2) + pow(rmu.y, 2) + pow(rmu.z, 2));

    rmu.x /= rm;
    rmu.y /= rm;
    rmu.z /= rm;

    rcu.x = coords[wi].x - topo.solvent_center.x;
    rcu.y = coords[wi].y - topo.solvent_center.y;
    rcu.z = coords[wi].z - topo.solvent_center.z;
    rc = sqrt(pow(rcu.x, 2) + pow(rcu.y, 2) + pow(rcu.z, 2));
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
    coord_t* coords, dvel_t* dvelocities, topo_t topo,
    double* theta, md_t md, double* energy,
    int* water_rank, int* water_shell) {
    int idx = blockIdx.x * blockDim.x + threadIdx.x;
    if (idx >= n_waters) return;

    int il = water_rank[idx];
    int is = water_shell[idx];
    if (is < 0) return;  // water not in any shell this step

    int wi, ii;
    coord_t rmu, rcu, f1O, f1H1, f1H2, f2;
    double rm, rc;
    double cos_th;
    double avtdum, arg, f0, dv;
    double ener;

    avtdum = 0;
    ii = idx;
    arg = 1 + ((1 - 2 * (double)(il + 1)) / (double)wshells[is].n_inshell);
    double theta_val = acos(arg);
    theta_val = theta_val - 3 * sin(theta_val) * wshells[is].cstb / 2;
    if (theta_val < 0) theta_val = 0;
    if (theta_val > M_PI) theta_val = M_PI;

    avtdum += theta[ii];
    ener = .5 * md.polarisation_force * pow(theta[ii] - theta_val + wshells[is].theta_corr, 2);
    // E_restraint.Upolx += ener;
    atomicAdd(energy, ener);

    dv = md.polarisation_force * (theta[ii] - theta_val + wshells[is].theta_corr);
    wi = n_atoms_solute + 3 * ii;

    rmu.x = coords[wi + 1].x + coords[wi + 2].x - 2 * coords[wi].x;
    rmu.y = coords[wi + 1].y + coords[wi + 2].y - 2 * coords[wi].y;
    rmu.z = coords[wi + 1].z + coords[wi + 2].z - 2 * coords[wi].z;

    rm = sqrt(pow(rmu.x, 2) + pow(rmu.y, 2) + pow(rmu.z, 2));

    rmu.x /= rm;
    rmu.y /= rm;
    rmu.z /= rm;

    rcu.x = coords[wi].x - topo.solvent_center.x;
    rcu.y = coords[wi].y - topo.solvent_center.y;
    rcu.z = coords[wi].z - topo.solvent_center.z;
    rc = sqrt(pow(rcu.x, 2) + pow(rcu.y, 2) + pow(rcu.z, 2));
    rcu.x /= rc;
    rcu.y /= rc;
    rcu.z /= rc;

    cos_th = rmu.x * rcu.x + rmu.y * rcu.y + rmu.z * rcu.z;
    if (cos_th > 1) cos_th = 1;
    if (cos_th < -1) cos_th = -1;
    f0 = sin(acos(cos_th));
    if (abs(f0) < 1.0E-12) f0 = 1.0E-12;
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

    atomicAdd(&dvelocities[wi].x, f0 * (f1O.x + f2.x));
    atomicAdd(&dvelocities[wi].y, f0 * (f1O.y + f2.y));
    atomicAdd(&dvelocities[wi].z, f0 * (f1O.z + f2.z));
    atomicAdd(&dvelocities[wi + 1].x, f0 * (f1H1.x));
    atomicAdd(&dvelocities[wi + 1].y, f0 * (f1H1.y));
    atomicAdd(&dvelocities[wi + 1].z, f0 * (f1H1.z));
    atomicAdd(&dvelocities[wi + 2].x, f0 * (f1H2.x));
    atomicAdd(&dvelocities[wi + 2].y, f0 * (f1H2.y));
    atomicAdd(&dvelocities[wi + 2].z, f0 * (f1H2.z));

    atomicAdd(&wshells[is].avtheta, avtdum / (double)wshells[is].n_inshell);
    atomicAdd(&wshells[is].avn_inshell, wshells[is].n_inshell);
}

void sort_waters() {
    using namespace CudaPolxWaterForce;

    int imin, jmin, jw;
    double tmin;
    // Sort the waters according to theta
    for (int is = 0; is < n_shells; is++) {
        imin = 0;
        for (int il = 0; il < wshells[is].n_inshell; il++) {
            tmin = 2 * M_PI;
            for (int jl = 0; jl < wshells[is].n_inshell; jl++) {
                // printf("Searching water %d in shell %d, total number: %d\n", jl, is, wshells[is].n_inshell);
                jw = polx_list_sh[jl * n_shells + is];
                // printf("zzzzzzzzz %d\n", jw);
                // printf("Water %d in shell %d has theta %f\n", jw, is, tdum[jw] * 180 / M_PI);
                if (tdum[jw] < tmin) {
                    jmin = jw;
                    tmin = tdum[jw];
                }
            }
            nsort[imin][is] = jmin;
            water_rank[jmin] = imin;
            water_shell[jmin] = is;
            imin++;
            tdum[jmin] = 99999;
        }
    }
}

void calc_polx_water_forces_host(int iteration) {
    CudaContext& ctx = CudaContext::instance();

    for (int is = 0; is < n_shells; is++) {
        wshells[is].n_inshell = 0;
        if (iteration == 0) {
            wshells[is].theta_corr = 0;
        }
    }
    using namespace CudaPolxWaterForce;
    if (!is_initialized) {
        water_rank = new int[n_waters];
        water_shell = new int[n_waters];
        polx_list_sh = new int[n_max_inshell * n_shells];

        is_initialized = true;
    }
    ctx.sync_all_to_device();

    coord_t* d_coords = ctx.d_coords;
    dvel_t* d_dvelocities = ctx.d_dvelocities;
    shell_t* d_wshells = ctx.d_wshells;
    int* d_list_sh = ctx.d_list_sh;
    double* d_theta = ctx.d_theta;
    double* d_theta0 = ctx.d_theta0;
    double* d_tdum = ctx.d_tdum;
    int* d_water_rank = ctx.d_water_rank;
    int* d_water_shell = ctx.d_water_shell;

    int blockSize = 256;
    int numBlocks = (n_waters + blockSize - 1) / blockSize;
    // printf("Calculated theta for %d waters in %d shells\n", n_waters, n_shells);
    calc_polx_theta_and_shells<<<numBlocks, blockSize>>>(
        n_waters, n_shells, n_atoms_solute, d_coords, topo, d_wshells, d_list_sh, d_theta, d_theta0, d_tdum);
    // printf("Calculated theta for %d waters in %d shells\n", n_waters, n_shells);

    // todo: sort in cpu now..
    cudaMemcpy(wshells, d_wshells, n_shells * sizeof(shell_t), cudaMemcpyDeviceToHost);
    cudaMemcpy(tdum, d_tdum, n_waters * sizeof(double), cudaMemcpyDeviceToHost);
    cudaMemcpy(polx_list_sh, d_list_sh, n_max_inshell * n_shells * sizeof(int), cudaMemcpyDeviceToHost);

    // Reset per-water metadata; only waters placed in shells will be overwritten in sort_waters().
    for (int i = 0; i < n_waters; i++) {
        water_rank[i] = -1;
        water_shell[i] = -1;
    }

    // printf("Sorting %d waters in %d shells\n", n_waters, n_shells);
    sort_waters();
    // printf("Sorted %d waters in %d shells\n", n_waters, n_shells);

    cudaMemcpy(d_water_rank, water_rank, n_waters * sizeof(int), cudaMemcpyHostToDevice);
    cudaMemcpy(d_water_shell, water_shell, n_waters * sizeof(int), cudaMemcpyHostToDevice);

    // Update theta_corr, averages
    if (iteration != 0 && iteration % itdis_update == 0) {
        for (int is = 0; is < n_shells; is++) {
            printf("SHELL %d\n", is);
            wshells[is].avtheta /= (double)itdis_update;
            wshells[is].avn_inshell /= (double)itdis_update;
            wshells[is].theta_corr = wshells[is].theta_corr + wshells[is].avtheta - acos(wshells[is].cstb);
            printf("average theta = %f, average in shell = %f, theta_corr = %f\n",
                   wshells[is].avtheta * 180 / M_PI, wshells[is].avn_inshell, wshells[is].theta_corr * 180 / M_PI);
            wshells[is].avtheta = 0;
            wshells[is].avn_inshell = 0;
        }
        cudaMemcpy(d_wshells, wshells, n_shells * sizeof(shell_t), cudaMemcpyHostToDevice);
    }

    // Calculate energy and force

    double* d_energy;
    check_cudaMalloc((void**)&d_energy, sizeof(double));
    calc_polx_water_forces_kernel<<<numBlocks, blockSize>>>(
        n_waters, n_atoms_solute, d_wshells, d_coords, d_dvelocities, topo,
        d_theta, md, d_energy, d_water_rank, d_water_shell);
    double energy;
    cudaMemcpy(&energy, d_energy, sizeof(double), cudaMemcpyDeviceToHost);
    E_restraint.Upolx += energy;
    cudaMemcpy(wshells, d_wshells, n_shells * sizeof(shell_t), cudaMemcpyDeviceToHost);
    // Copy back forces for all atoms (solute + solvent); water forces were being dropped.
    cudaMemcpy(dvelocities, d_dvelocities, n_atoms * sizeof(dvel_t), cudaMemcpyDeviceToHost);
    cudaFree(d_energy);
}

void cleanup_polx_water_force(

) {
    using namespace CudaPolxWaterForce;
    if (is_initialized) {
        delete[] water_rank;
        delete[] water_shell;

        is_initialized = false;
    }
}
