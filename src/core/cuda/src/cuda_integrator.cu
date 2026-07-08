#include "cuda_force_accumulation.cuh"
#include "cuda_integrator.cuh"
#include "temperature.h"

namespace {
__global__ void calc_leapfrog_kernel(
    atype_t* atypes,
    catype_t* catypes,
    vel_t* velocities,
    dvel_t* dvelocities,
    coord_t* coords,
    coord_t* xcoords,
    int n_atoms,
    int n_atoms_solute,
    const double* temperature_results,
    double dt) {
    int i = blockIdx.x * blockDim.x + threadIdx.x;
    if (i >= n_atoms) return;

    const double mass_i = catypes[atypes[i].code - 1].m;
    const double winv_i = 1 / mass_i;
    const double scale = (i < n_atoms_solute) ? temperature_results[R_TSCALE_SOL]
                                              : temperature_results[R_TSCALE_SLV];

    const double fx = force_from_accum(dvelocities[i].x);
    const double fy = force_from_accum(dvelocities[i].y);
    const double fz = force_from_accum(dvelocities[i].z);

    velocities[i].x = (velocities[i].x - fx * dt * winv_i) * scale;
    velocities[i].y = (velocities[i].y - fy * dt * winv_i) * scale;
    velocities[i].z = (velocities[i].z - fz * dt * winv_i) * scale;

    xcoords[i].x = coords[i].x;
    xcoords[i].y = coords[i].y;
    xcoords[i].z = coords[i].z;

    coords[i].x += velocities[i].x * dt;
    coords[i].y += velocities[i].y * dt;
    coords[i].z += velocities[i].z * dt;
}

__global__ void update_velocities_from_positions_kernel(
    vel_t* velocities,
    const coord_t* coords,
    const coord_t* xcoords,
    int n_atoms,
    double dt) {
    int i = blockIdx.x * blockDim.x + threadIdx.x;
    if (i >= n_atoms) return;

    velocities[i].x = (coords[i].x - xcoords[i].x) / dt;
    velocities[i].y = (coords[i].y - xcoords[i].y) / dt;
    velocities[i].z = (coords[i].z - xcoords[i].z) / dt;
}

}  // namespace

void CudaIntegrator::step(Context& ctx) {
    auto d_atypes = ctx.atypes->gpu_data_p;
    auto d_catypes = ctx.catypes->gpu_data_p;
    auto d_velocities = ctx.velocities->gpu_data_p;
    auto d_dvelocities = ctx.dvelocities->gpu_data_p;
    auto d_coords = ctx.coords->gpu_data_p;
    auto d_xcoords = ctx.xcoords->gpu_data_p;

    const double* d_temperature_results = temperature_->data().results->gpu_data_p;

    int blockSize = 256;
    int numBlocks = (ctx.n_atoms + blockSize - 1) / blockSize;
    calc_leapfrog_kernel<<<numBlocks, blockSize>>>(
        d_atypes, d_catypes, d_velocities, d_dvelocities,
        d_coords, d_xcoords, ctx.n_atoms, ctx.n_atoms_solute,
        d_temperature_results, ctx.dt);

    if (shake_->data().n_constraints > 0) {
        shake_->apply(ctx);
        update_velocities_from_positions_kernel<<<numBlocks, blockSize>>>(
            d_velocities, d_coords, d_xcoords, ctx.n_atoms, ctx.dt);
    }
}