#include "cuda_nonbonded_force.cuh"
#include "cuda_nonbonded_pair.cuh"

namespace {

__global__ void update_nonbonded_coords_kernel(
    const coord_t* coords, const int* atom_idx,
    real_t* cx, real_t* cy, real_t* cz, int sz) {
    const int i = blockIdx.x * blockDim.x + threadIdx.x;
    if (i >= sz) return;
    const int idx = atom_idx[i];
    if (idx < 0) {  // padding slot (atom_idx == -1); main kernel treats these as empty
        cx[i] = cy[i] = cz[i] = 0;
        return;
    }
    cx[i] = static_cast<real_t>(coords[idx].x);
    cy[i] = static_cast<real_t>(coords[idx].y);
    cz[i] = static_cast<real_t>(coords[idx].z);
}

__global__ void nonbonded_kernel(
    // ---- dimensions ----
    int sz,        // data_.n_total, number of participating atoms
    int n_states,  // ctx.n_lambdas, used by nb_coul_slot

    // ---- per-atom arrays (length sz, parallel to atom_idx) ----
    const int* atom_idx,               // data_.atom_idx, local i -> global atom index
    const uint8_t* category,           // data_.category, P/Q/W
    const int* q_state,                // data_.q_state, Q state; -1 for P/W
    const real_t* atom_lambdas,        // data_.atom_lambdas
    const real_t* atom_charge,         // data_.atom_charge
    const vdw_atom_param_t* atom_vdw,  // data_.atom_vdw

    // ---- exclusion data ----
    const uint32_t* pair_masks,

    // ---- topology scalars (passed by value) ----
    real_t coulomb_constant,  // ctx.topo.coulomb_constant
    int vdw_rule,

    // ---- coordinates / outputs ----
    const real_t* cx, const real_t* cy, const real_t* cz,
    dvel_t* dvelocities,  // ctx.dvelocities->gpu_data_p (fixed-point, atomic_add_force)

    // ---- energy accumulators  ----
    energy_accum_t* e) {
    const int block_num = (sz + 31) >> 5;
    const int total_tiles = (block_num * (block_num + 1)) >> 1;
    const int warps_per_block = blockDim.x >> 5;
    const int tid = threadIdx.x;
    const int lane = tid & 31;
    const int warp_in_block = tid >> 5;

    const int tile = blockIdx.x * warps_per_block + warp_in_block;
    if (tile >= total_tiles) return;
    const size_t row = static_cast<size_t>(tile) * 32 + lane;
    const uint32_t active_mask = pair_masks[row];

    auto [tile_x, tile_y] = get_tile_idx(block_num, tile);

    const int base_x = tile_x << 5;
    const int base_y = tile_y << 5;

    int x_idx = base_x + lane;
    int y_idx = base_y + lane;
    x_idx = x_idx < sz ? x_idx : -1;
    y_idx = y_idx < sz ? y_idx : -1;

    nonbonded_force_calculation(x_idx, y_idx, active_mask, base_x, base_y, n_states,
                                atom_idx, category, q_state, atom_lambdas, atom_charge, atom_vdw, coulomb_constant, vdw_rule, cx, cy, cz, dvelocities, e);
}

}  // namespace

void CudaNonbondedForce::init_backend(Context& ctx) {
    // Buffers are indexed by combined-list position [0, n_total), which exceeds
    // n_atoms because Q atoms are duplicated per FEP state and the list is padded.
    coord_x_ = std::make_unique<HostDeviceBuffer<real_t>>(data_.n_total);
    coord_y_ = std::make_unique<HostDeviceBuffer<real_t>>(data_.n_total);
    coord_z_ = std::make_unique<HostDeviceBuffer<real_t>>(data_.n_total);

    lrf_.reset();
    if (ctx.md.lrf && ctx.md.non_bond != 0) {
        lrf_ = std::make_unique<CudaLrfForce>();
        lrf_->init(ctx, data_);
    }
}

void CudaNonbondedForce::calc_all_direct_pairs(Context& ctx) {
    const int thread_num = 256;
    int tile_num_per_block = thread_num >> 5;
    int n_atom = data_.n_total;
    int block_num = (n_atom + 31) >> 5;
    int total_tiles = block_num * (block_num + 1) >> 1;
    int grid_sz = (total_tiles + tile_num_per_block - 1) / tile_num_per_block;

    dim3 grid = dim3(grid_sz);
    nonbonded_kernel<<<grid, thread_num>>>(n_atom, ctx.n_lambdas(),
                                           data_.atom_idx->gpu_data_p, data_.category->gpu_data_p, data_.q_state->gpu_data_p,
                                           data_.atom_lambdas->gpu_data_p, data_.atom_charge->gpu_data_p, data_.atom_vdw->gpu_data_p,
                                           data_.pair_masks->gpu_data_p, ctx.topo.coulomb_constant, ctx.topo.vdw_rule,
                                           coord_x_->gpu_data_p, coord_y_->gpu_data_p, coord_z_->gpu_data_p, ctx.dvelocities->gpu_data_p, ctx.energy.device());
}

void CudaNonbondedForce::calc(Context& ctx) {
    const int size = data_.n_total;
    if (size == 0) return;

    constexpr int threads = 256;

    update_nonbonded_coords_kernel<<<(size + threads - 1) / threads, threads>>>(
        ctx.coords->gpu_data_p,
        data_.atom_idx->gpu_data_p,
        coord_x_->gpu_data_p,
        coord_y_->gpu_data_p,
        coord_z_->gpu_data_p,
        size);

    check_cuda(cudaGetLastError());

    if (lrf_) {
        lrf_->calc(ctx, {coord_x_->gpu_data_p, coord_y_->gpu_data_p, coord_z_->gpu_data_p});
    } else {
        calc_all_direct_pairs(ctx);
    }
}