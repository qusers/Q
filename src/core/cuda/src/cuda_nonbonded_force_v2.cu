#include "cuda_nonbonded_force_v2.cuh"

namespace {
__device__ int2 get_tile_idx(int n, int t) {
    int x = (int)floorf((2 * n + 1 - sqrtf((2 * n + 1) * (2 * n + 1) - 8 * t)) * 0.5f);
    int y = t - (x * n - (x * (x - 1) >> 1));
    if (y < 0) {
        x--;
        y += (n - x);
    }
    y += x;
    return {x, y};
}

}  // namespace

void CudaNonbondedForce::init_backend(Context& ctx) {
    int n_slots = nb_total_slots(ctx.n_lambdas);
    e_coul_ = std::make_unique<HostDeviceBuffer<energy_accum_t>>(n_slots);
    e_vdw_ = std::make_unique<HostDeviceBuffer<energy_accum_t>>(n_slots);
}

__global__ void nonbonded_kernel(
    // ---- dimensions ----
    int sz,              // data_.n_total, number of participating atoms
    int n_states,        // ctx.n_lambdas, used by nb_energy_slot
    int n_charge_types,  // data_.n_charge_types, row stride of charge table
    int n_catype_types,  // data_.n_catype_types, row stride of catype table
    int n_atoms_solute,  // ctx.n_atoms_solute, water grouping + LJ_matrix row stride

    // ---- per-atom arrays (length sz, parallel to atom_idx) ----
    const int* atom_idx,         // data_.atom_idx, local i -> global atom index
    const uint8_t* category,     // data_.category, P/Q/W
    const int* q_state,          // data_.q_state, Q state; -1 for P/W
    const real_t* atom_lambdas,  // data_.atom_lambdas
    const int* charge_types,     // data_.charge_types
    const int* catype_types,     // data_.catype_types

    // ---- precomputed pairwise tables ----
    const real_t* charge_pair_products,          // data_.charge_pair_products [n_ct*n_ct]
    const vdw_pair_param_t* catype_pair_params,  // data_.catype_pair_params  [n_cat*n_cat]

    // ---- exclusion data ----
    const int* LJ_matrix,  // ctx.LJ_matrix->gpu_data_p

    // ---- topology scalars (passed by value) ----
    real_t el14_scale,        // ctx.topo.el14_scale
    real_t coulomb_constant,  // ctx.topo.coulomb_constant

    // ---- coordinates / outputs ----
    const nonbond_coord_t* coords,  // ctx.nonbond_coord_t->gpu_data_p
    dvel_t* dvelocities,            // ctx.dvelocities->gpu_data_p (fixed-point, atomic_add_force)

    // ---- energy accumulators (length nb_total_slots(n_states)) ----
    energy_accum_t* e_coul,  // e_coul_->gpu_data_p
    energy_accum_t* e_vdw    // e_vdw_->gpu_data_p
) {
    const int block_num = (sz + 31) >> 5;
    const int total_tiles = (block_num * (block_num + 1)) >> 1;
    const int warps_per_block = blockDim.x >> 5;
    const int tid = threadIdx.x;
    const int lane = tid & 31;
    const int warp_in_block = tid >> 5;

    const int tile = blockIdx.x * warps_per_block + warp_in_block;
    if (tile >= total_tiles) return;

    auto [tile_x, tile_y] = get_tile_idx(block_num, tile);

    const int base_x = tile_x << 5;
    const int base_y = tile_y << 5;

    int x_idx = base_x + lane;
    int y_idx = base_y + lane;
}

void CudaNonbondedForce::calc(Context& ctx) {
    const int thread_num = 256;

    int tile_num_per_block = thread_num >> 5;
    int n_atom = data_.n_total;
    int block_num = (n_atom + 31) >> 5;
    int total_tiles = block_num * (block_num + 1) >> 1;
    int grid_sz = (total_tiles + tile_num_per_block - 1) / tile_num_per_block;

    auto& d_e_coul = e_coul_->gpu_data_p;
    auto& d_e_vdw = e_vdw_->gpu_data_p;

    dim3 grid = dim3(grid_sz);
    cudaMemset(d_e_coul, 0, sizeof(energy_accum_t) * e_coul_->length);
    cudaMemset(d_e_vdw, 0, sizeof(energy_accum_t) * e_vdw_->length);
}