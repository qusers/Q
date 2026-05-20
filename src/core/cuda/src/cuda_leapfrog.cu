#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <iostream>

#include "common/include/context.h"
#include "cuda/include/cuda_force_accum.cuh"
#include "cuda/include/cuda_leapfrog.cuh"
#include "cuda/include/cuda_shake_constraints.cuh"
#include "cuda_utility.cuh"

namespace CudaLeapfrog {

int debug_leapfrog_step = 0;

real_t debug_force_component(force_fixed_storage_t value) {
    return fixed_to_force(force_fixed_from_storage(value));
}

}

__global__ void calc_leapfrog_kernel(
    atype_t* atypes,
    catype_t* catypes,
    vel_t* velocities,
    cuda_dvel_t* dvelocities,
    coord_t* coords,
    coord_t* xcoords,
    int n_atoms,
    int n_atoms_solute,
    real_t Tscale_solute,
    real_t Tscale_solvent,
    real_t dt) {
    int idx = blockIdx.x * blockDim.x + threadIdx.x;
    if (idx >= n_atoms) return;
    int i = idx;

    // Kernel implementation goes here
    real_t mass_i, winv_i;

    mass_i = catypes[atypes[i].code - 1].m;

    winv_i = 1 / mass_i;
    real_t scale = (i < n_atoms_solute) ? Tscale_solute : Tscale_solvent;
    const real_t dvel_x = force_accum_component_value(dvelocities[i].x);
    const real_t dvel_y = force_accum_component_value(dvelocities[i].y);
    const real_t dvel_z = force_accum_component_value(dvelocities[i].z);
    velocities[i].x = (velocities[i].x - dvel_x * dt * winv_i) * scale;
    velocities[i].y = (velocities[i].y - dvel_y * dt * winv_i) * scale;
    velocities[i].z = (velocities[i].z - dvel_z * dt * winv_i) * scale;

    xcoords[i].x = coords[i].x;
    xcoords[i].y = coords[i].y;
    xcoords[i].z = coords[i].z;

    coords[i].x += velocities[i].x * dt;
    coords[i].y += velocities[i].y * dt;
    coords[i].z += velocities[i].z * dt;
}

void calc_leapfrog_host() {
    auto& host = Context::instance();
    auto d_atypes = host.atypes->gpu_data_p;
    auto d_catypes = host.catypes->gpu_data_p;
    auto d_velocities = host.velocities->gpu_data_p;
    auto d_dvelocities = cuda_force_accum_buffer(host);
    auto d_coords = host.coords->gpu_data_p;
    auto d_xcoords = host.xcoords->gpu_data_p;

    int blockSize = 256;
    int numBlocks = (host.n_atoms + blockSize - 1) / blockSize;
    calc_leapfrog_kernel<<<numBlocks, blockSize>>>(
        d_atypes,
        d_catypes,
        d_velocities,
        d_dvelocities,
        d_coords,
        d_xcoords,
        host.n_atoms,
        host.n_atoms_solute,
        host.Tscale_solute,
        host.Tscale_solvent,
        host.dt);
    check_cuda(cudaDeviceSynchronize());

    host.velocities->download();
    host.coords->download();
    host.xcoords->download();

    CudaLeapfrog::debug_leapfrog_step++;
    if (const char* atoms_env = std::getenv("QGPU_DEBUG_SHAKE_ATOMS")) {
        int debug_atoms[3] = {0, 0, 0};
        if (std::sscanf(atoms_env, "%d,%d,%d", &debug_atoms[0], &debug_atoms[1], &debug_atoms[2]) == 3) {
            int debug_start = 0;
            if (const char* start_env = std::getenv("QGPU_DEBUG_SHAKE_START")) {
                debug_start = std::atoi(start_env);
            }
            int valid = 1;
            for (int k = 0; k < 3; k++) {
                debug_atoms[k]--;
                if (debug_atoms[k] < 0 || debug_atoms[k] >= host.n_atoms) valid = 0;
            }
            if (valid) {
                cuda_dvel_t debug_dvel[3] = {};
                check_cuda(cudaMemcpy(&debug_dvel[0], d_dvelocities + debug_atoms[0], sizeof(cuda_dvel_t), cudaMemcpyDeviceToHost));
                check_cuda(cudaMemcpy(&debug_dvel[1], d_dvelocities + debug_atoms[1], sizeof(cuda_dvel_t), cudaMemcpyDeviceToHost));
                check_cuda(cudaMemcpy(&debug_dvel[2], d_dvelocities + debug_atoms[2], sizeof(cuda_dvel_t), cudaMemcpyDeviceToHost));

                auto* coords = host.coords->cpu_data_p;
                auto* xcoords = host.xcoords->cpu_data_p;
                auto* velocities = host.velocities->cpu_data_p;
                const int pairs[3][2] = {{0, 1}, {0, 2}, {1, 2}};
                bool print_debug = CudaLeapfrog::debug_leapfrog_step >= debug_start;
                double dist[3] = {};
                double old_dist[3] = {};
                double scp[3] = {};
                for (int p = 0; p < 3; p++) {
                    const int ai = debug_atoms[pairs[p][0]];
                    const int aj = debug_atoms[pairs[p][1]];
                    const double dx = static_cast<double>(coords[ai].x) - static_cast<double>(coords[aj].x);
                    const double dy = static_cast<double>(coords[ai].y) - static_cast<double>(coords[aj].y);
                    const double dz = static_cast<double>(coords[ai].z) - static_cast<double>(coords[aj].z);
                    const double odx = static_cast<double>(xcoords[ai].x) - static_cast<double>(xcoords[aj].x);
                    const double ody = static_cast<double>(xcoords[ai].y) - static_cast<double>(xcoords[aj].y);
                    const double odz = static_cast<double>(xcoords[ai].z) - static_cast<double>(xcoords[aj].z);
                    dist[p] = std::sqrt(dx * dx + dy * dy + dz * dz);
                    old_dist[p] = std::sqrt(odx * odx + ody * ody + odz * odz);
                    scp[p] = dx * odx + dy * ody + dz * odz;
                    if (dist[p] > 2.0 || std::fabs(scp[p]) < 0.05) print_debug = true;
                }
                if (print_debug) {
                    std::printf(">>> DEBUG pre-shake step=%d atoms=%d,%d,%d\n",
                                CudaLeapfrog::debug_leapfrog_step,
                                debug_atoms[0] + 1, debug_atoms[1] + 1, debug_atoms[2] + 1);
                    for (int k = 0; k < 3; k++) {
                        const int ai = debug_atoms[k];
                        std::printf(">>> DEBUG atom %d coord=(%.9g %.9g %.9g) xcoord=(%.9g %.9g %.9g) vel=(%.9g %.9g %.9g) dvel=(%.9g %.9g %.9g)\n",
                                    ai + 1,
                                    static_cast<double>(coords[ai].x), static_cast<double>(coords[ai].y), static_cast<double>(coords[ai].z),
                                    static_cast<double>(xcoords[ai].x), static_cast<double>(xcoords[ai].y), static_cast<double>(xcoords[ai].z),
                                    static_cast<double>(velocities[ai].x), static_cast<double>(velocities[ai].y), static_cast<double>(velocities[ai].z),
                                    static_cast<double>(CudaLeapfrog::debug_force_component(debug_dvel[k].x)),
                                    static_cast<double>(CudaLeapfrog::debug_force_component(debug_dvel[k].y)),
                                    static_cast<double>(CudaLeapfrog::debug_force_component(debug_dvel[k].z)));
                    }
                    for (int p = 0; p < 3; p++) {
                        std::printf(">>> DEBUG pair %d-%d dist=%.9g old_dist=%.9g scp=%.9g\n",
                                    debug_atoms[pairs[p][0]] + 1, debug_atoms[pairs[p][1]] + 1, dist[p], old_dist[p], scp[p]);
                    }
                }
            }
        }
    }

    // shake
    printf("n_shake_constraints: %d\n", host.n_shake_constraints);
    if (host.n_shake_constraints > 0) {
        calc_shake_constraints_host();
        auto& velocities = host.velocities->cpu_data_p;
        auto& coords = host.coords->cpu_data_p;
        auto* xcoords = host.xcoords->cpu_data_p;
        for (int i = 0; i < host.n_atoms; i++) {
            velocities[i].x = (coords[i].x - xcoords[i].x) / host.dt;
            velocities[i].y = (coords[i].y - xcoords[i].y) / host.dt;
            velocities[i].z = (coords[i].z - xcoords[i].z) / host.dt;
        }
        host.velocities->upload();
    }
}

void init_leapfrog_kernel_data() {}

void cleanup_leapfrog() {}
