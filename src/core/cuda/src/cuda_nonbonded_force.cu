#include "cuda_force_accumulation.cuh"
#include "cuda_nonbonded_force.cuh"

namespace {
// In the continuous 16-lane, left rotate one
// In the current lane, read (lane + 1) % 16
// All lane group should be active
#if defined(QGPU_BACKEND_HIP)
template <typename T>
__device__ __forceinline__ T rotate_left16_dpp(T value) {
    static_assert(sizeof(T) == 4 || sizeof(T) == 8, "Only 32-bit and 64-bit values are supported");

    if constexpr (sizeof(T) == 4) {
        int bits = __builtin_bit_cast(int, value);

        bits = __builtin_amdgcn_mov_dpp( bits, 0x12f, 0xf, 0xf, false);

        return __builtin_bit_cast(T, bits);
    } else {
        struct Words {
            int word0;
            int word1;
        };
        static_assert(sizeof(Words) == sizeof(T));

        Words bits = __builtin_bit_cast(Words, value);

        bits.word0 = __builtin_amdgcn_mov_dpp(bits.word0, 0x12f, 0xf, 0xf, false);
        bits.word1 = __builtin_amdgcn_mov_dpp(bits.word1, 0x12f, 0xf, 0xf, false);

        return __builtin_bit_cast(T, bits);
    }
}

// In the 32-lane group, switch the left half group and the right half group 
// 0 <-> 16, 1 <-> 17, ..., 15 <-> 31。
template <typename T>
__device__ __forceinline__
T swap_halves32_value(T value) {
    return __shfl_xor(value, 16, 32);
}


__device__ __forceinline__ void shuffle_left_16_dpp(
    int& atom,
    real_t& atom_charge,
    vdw_atom_param_t& atom_vdw,
    real_t3& atom_force,
    real_t3& atom_coord)
{
    atom = rotate_left16_dpp(atom);
    atom_charge = rotate_left16_dpp(atom_charge);

    atom_vdw.aii_normal = rotate_left16_dpp(atom_vdw.aii_normal);
    atom_vdw.bii_normal = rotate_left16_dpp(atom_vdw.bii_normal);
    atom_vdw.aii_14 = rotate_left16_dpp(atom_vdw.aii_14);
    atom_vdw.bii_14 = rotate_left16_dpp(atom_vdw.bii_14);

    atom_force.x = rotate_left16_dpp(atom_force.x);
    atom_force.y = rotate_left16_dpp(atom_force.y);
    atom_force.z = rotate_left16_dpp(atom_force.z);

    atom_coord.x = rotate_left16_dpp(atom_coord.x);
    atom_coord.y = rotate_left16_dpp(atom_coord.y);
    atom_coord.z = rotate_left16_dpp(atom_coord.z);
}

__device__ __forceinline__
void swap_halves_32(
    int& atom,
    real_t& atom_charge,
    vdw_atom_param_t& atom_vdw,
    real_t3& atom_force,
    real_t3& atom_coord)
{
    atom = swap_halves32_value(atom);
    atom_charge = swap_halves32_value(atom_charge);

    atom_vdw.aii_normal = swap_halves32_value(atom_vdw.aii_normal);
    atom_vdw.bii_normal = swap_halves32_value(atom_vdw.bii_normal);
    atom_vdw.aii_14 = swap_halves32_value(atom_vdw.aii_14);
    atom_vdw.bii_14 = swap_halves32_value(atom_vdw.bii_14);

    atom_force.x = swap_halves32_value(atom_force.x);
    atom_force.y = swap_halves32_value(atom_force.y);
    atom_force.z = swap_halves32_value(atom_force.z);

    atom_coord.x = swap_halves32_value(atom_coord.x);
    atom_coord.y = swap_halves32_value(atom_coord.y);
    atom_coord.z = swap_halves32_value(atom_coord.z);
}
#endif

template <typename T>
__device__ __forceinline__
T shfl32(T value, int src_lane) {
#if defined(QGPU_BACKEND_HIP)
    return __shfl(value, src_lane, 32);
#else
    return __shfl_sync(0xffffffffu, value, src_lane, 32);
#endif
}

template <typename T>
__device__ __forceinline__
T shfl_down32(T value, unsigned delta) {
#if defined(QGPU_BACKEND_HIP)
    return __shfl_down(value, delta, 32);
#else
    return __shfl_down_sync(0xffffffffu, value, delta, 32);
#endif
}

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

__device__ void compute_pair(
    // atom1
    int atom1, uint8_t atom1_type, int atom1_state,
    real_t atom1_charge, const vdw_atom_param_t& atom1_vdw, real_t atom1_lambda, const real_t3& atom1_coord,
    /// atom2
    int atom2, uint8_t atom2_type, int atom2_state,
    real_t atom2_charge, const vdw_atom_param_t& atom2_vdw, real_t atom2_lambda, const real_t3& atom2_coord,
    // exclusion
    int n_atoms_solute, const int* LJ_matrix,
    // scalars
    real_t el14_scale, real_t coulomb_constant, int vdw_rule,
    // output
    real_t3& atom1_force, real_t3& atom2_force,
    real_t& e_coul, real_t& e_vdw) {
    if (atom1 == -1 || atom2 == -1 || atom1 == atom2) return;
    auto bond_type = get_bond_type(n_atoms_solute, LJ_matrix, atom1, atom1_type, atom2, atom2_type);
    if (bond_type == BondType::Bond23) return;

    constexpr uint8_t Q = static_cast<uint8_t>(AtomCategory::Q);
    if (atom1_type == Q && atom2_type == Q && atom1_state != atom2_state) return;

    real_t dx = atom2_coord.x - atom1_coord.x;
    real_t dy = atom2_coord.y - atom1_coord.y;
    real_t dz = atom2_coord.z - atom1_coord.z;
    real_t dis2 = dx * dx + dy * dy + dz * dz;
    real_t inv_dis = rsqrt(dis2);

    bool is_14 = (bond_type == BondType::Bond14);
    real_t qij = atom1_charge * atom2_charge;
    real_t scaling = is_14 ? el14_scale : 1;
    real_t2 pair = (is_14) ? combine_vdw(vdw_rule, atom1_vdw.aii_14, atom1_vdw.bii_14, atom2_vdw.aii_14, atom2_vdw.bii_14) : combine_vdw(vdw_rule, atom1_vdw.aii_normal, atom1_vdw.bii_normal, atom2_vdw.aii_normal, atom2_vdw.bii_normal);

    auto [vel, dvel] = calc_electrostatic(qij * scaling, coulomb_constant, inv_dis);
    auto [vvdw, dvvdw] = calc_vdw(pair, inv_dis);

    real_t lambda = min(atom1_lambda, atom2_lambda);

    real_t dva = (dvel + dvvdw) * inv_dis * lambda;

    atom1_force.x -= dva * dx;
    atom1_force.y -= dva * dy;
    atom1_force.z -= dva * dz;

    atom2_force.x += dva * dx;
    atom2_force.y += dva * dy;
    atom2_force.z += dva * dz;

    e_coul += vel;
    e_vdw += vvdw;
}

__device__ void shuffle(int& atom, real_t& atom_charge, vdw_atom_param_t& atom_vdw, real_t3& atom_force, real_t3& atom_coord) {
    int src = ((threadIdx.x & 31) + 1) & 31;
    atom = shfl32(atom, src);
    atom_charge = shfl32(atom_charge, src);

    atom_vdw.aii_normal = shfl32(atom_vdw.aii_normal, src);
    atom_vdw.bii_normal = shfl32(atom_vdw.bii_normal, src);
    atom_vdw.aii_14 = shfl32(atom_vdw.aii_14, src);
    atom_vdw.bii_14 = shfl32(atom_vdw.bii_14, src);

    atom_force.x = shfl32(atom_force.x, src);
    atom_force.y = shfl32(atom_force.y, src);
    atom_force.z = shfl32(atom_force.z, src);
    atom_coord.x = shfl32(atom_coord.x, src);
    atom_coord.y = shfl32(atom_coord.y, src);
    atom_coord.z = shfl32(atom_coord.z, src);
}

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

}  // namespace

void CudaNonbondedForce::init_backend(Context& ctx) {
    // Buffers are indexed by combined-list position [0, n_total), which exceeds
    // n_atoms because Q atoms are duplicated per FEP state and the list is padded.
    coord_x = std::make_unique<HostDeviceBuffer<real_t>>(data_.n_total);
    coord_y = std::make_unique<HostDeviceBuffer<real_t>>(data_.n_total);
    coord_z = std::make_unique<HostDeviceBuffer<real_t>>(data_.n_total);
}

__global__ void nonbonded_kernel(
    // ---- dimensions ----
    int sz,              // data_.n_total, number of participating atoms
    int n_states,        // ctx.n_lambdas, used by nb_coul_slot
    int n_atoms_solute,  // ctx.n_atoms_solute, water grouping + LJ_matrix row stride

    // ---- per-atom arrays (length sz, parallel to atom_idx) ----
    const int* atom_idx,               // data_.atom_idx, local i -> global atom index
    const uint8_t* category,           // data_.category, P/Q/W
    const int* q_state,                // data_.q_state, Q state; -1 for P/W
    const real_t* atom_lambdas,        // data_.atom_lambdas
    const real_t* atom_charge,         // data_.atom_charge
    const vdw_atom_param_t* atom_vdw,  // data_.atom_vdw

    // ---- exclusion data ----
    const int* LJ_matrix,  // ctx.LJ_matrix->gpu_data_p

    // ---- topology scalars (passed by value) ----
    real_t el14_scale,        // ctx.topo.el14_scale
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

    const auto tile_idx = get_tile_idx(block_num, tile);
    const int tile_x = tile_idx.x;
    const int tile_y = tile_idx.y;

    const int base_x = tile_x << 5;
    const int base_y = tile_y << 5;

    int x_idx = base_x + lane;
    int y_idx = base_y + lane;

    const int atom1 = x_idx < sz ? atom_idx[x_idx] : -1;

    const auto& atom1_type = atom1 == -1 ? static_cast<uint8_t>(AtomCategory::INVALID) : category[x_idx];
    const int atom1_state = atom1 == -1 ? -1 : q_state[x_idx];
    const real_t atom1_charge = atom1 == -1 ? 0 : atom_charge[x_idx];
    const vdw_atom_param_t atom1_vdw = atom1 == -1 ? vdw_atom_param_t{0, 0, 0, 0} : atom_vdw[x_idx];
    const real_t atom1_lambda = atom1 == -1 ? 0 : atom_lambdas[x_idx];
    real_t3 atom1_coord = atom1 == -1 ? real_t3{0, 0, 0} : real_t3{cx[x_idx], cy[x_idx], cz[x_idx]};
    real_t3 atom1_force = {0, 0, 0};

    int atom2 = y_idx < sz ? atom_idx[y_idx] : -1;
    const uint8_t tile_y_type = category[base_y];
    const int tile_y_state = q_state[base_y];
    const real_t tile_y_lambda = atom_lambdas[base_y];

    real_t atom2_charge = atom2 == -1 ? 0 : atom_charge[y_idx];
    vdw_atom_param_t atom2_vdw = atom2 == -1 ? vdw_atom_param_t{0, 0, 0, 0} : atom_vdw[y_idx];
    real_t3 atom2_coord = atom2 == -1 ? real_t3{0, 0, 0} : real_t3{cx[y_idx], cy[y_idx], cz[y_idx]};
    real_t3 atom2_force = {0, 0, 0};

    real_t local_e_coul = 0, local_e_vdw = 0;
    bool is_diag = (tile_x == tile_y);
    for (int i = 0; i < 32; i++) {
        if (!is_diag || atom1 < atom2) {
            compute_pair(atom1, atom1_type, atom1_state, atom1_charge, atom1_vdw, atom1_lambda, atom1_coord,
                         atom2, tile_y_type, tile_y_state, atom2_charge, atom2_vdw, tile_y_lambda, atom2_coord,
                         n_atoms_solute, LJ_matrix,
                         el14_scale, coulomb_constant, vdw_rule,
                         atom1_force, atom2_force,
                         local_e_coul, local_e_vdw);
        }
        #if defined(QGPU_BACKEND_HIP)
        shuffle_left_16_dpp( atom2, atom2_charge, atom2_vdw, atom2_force, atom2_coord);

        // 第 16、32 次计算后交换半组。
        // 第一次进入另一半组；第二次恢复原始分布。
        if ((i & 15) == 15) {
            swap_halves_32( atom2, atom2_charge, atom2_vdw, atom2_force, atom2_coord);
        }
        #else
            shuffle( atom2, atom2_charge, atom2_vdw, atom2_force, atom2_coord);
        #endif        


    }

    if (atom1 >= 0) {
        atomic_add_force(&dvelocities[atom1].x, atom1_force.x);
        atomic_add_force(&dvelocities[atom1].y, atom1_force.y);
        atomic_add_force(&dvelocities[atom1].z, atom1_force.z);
    }

    if (atom2 >= 0) {
        atomic_add_force(&dvelocities[atom2].x, atom2_force.x);
        atomic_add_force(&dvelocities[atom2].y, atom2_force.y);
        atomic_add_force(&dvelocities[atom2].z, atom2_force.z);
    }
    for (int offset = 16; offset > 0; offset >>= 1) {
        local_e_coul += shfl_down32(local_e_coul, offset);
        local_e_vdw += shfl_down32(local_e_vdw, offset);
    }

    if (lane == 0) {
        uint8_t tile_cat_x = category[base_x];
        uint8_t tile_cat_y = category[base_y];
        int tile_state_x = q_state[base_x];
        int tile_state_y = q_state[base_y];
        int coul_slot = nb_coul_slot(tile_cat_x, tile_cat_y, tile_state_x, tile_state_y, n_states);

        atomic_add_energy(&e[coul_slot], local_e_coul);
        atomic_add_energy(&e[coul_slot + 1], local_e_vdw);
    }
}

void CudaNonbondedForce::calc(Context& ctx) {
    /*
    Sync the coords to CudaNonbondedForce::coords first.
    */
    int sz = data_.n_total;
    int sync_block = 256;
    int sync_grid = (sz + sync_block - 1) / sync_block;
    update_nonbonded_coords_kernel<<<sync_grid, sync_block>>>(ctx.coords->gpu_data_p, data_.atom_idx->gpu_data_p, coord_x->gpu_data_p, coord_y->gpu_data_p, coord_z->gpu_data_p, sz);

    /*
    Do calculation
    */
    const int thread_num = 256;
    int tile_num_per_block = thread_num >> 5;
    int n_atom = data_.n_total;
    int block_num = (n_atom + 31) >> 5;
    int total_tiles = block_num * (block_num + 1) >> 1;
    int grid_sz = (total_tiles + tile_num_per_block - 1) / tile_num_per_block;

    dim3 grid = dim3(grid_sz);
    nonbonded_kernel<<<grid, thread_num>>>(n_atom, ctx.n_lambdas(), ctx.n_atoms_solute,
                                           data_.atom_idx->gpu_data_p, data_.category->gpu_data_p, data_.q_state->gpu_data_p,
                                           data_.atom_lambdas->gpu_data_p, data_.atom_charge->gpu_data_p, data_.atom_vdw->gpu_data_p,
                                           ctx.LJ_matrix->gpu_data_p, ctx.topo.el14_scale, ctx.topo.coulomb_constant, ctx.topo.vdw_rule,
                                           coord_x->gpu_data_p, coord_y->gpu_data_p, coord_z->gpu_data_p, ctx.dvelocities->gpu_data_p, ctx.energy.device());
}