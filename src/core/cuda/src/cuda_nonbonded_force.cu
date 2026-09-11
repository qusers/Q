#include <cub/device/device_scan.cuh>

#include "cuda_force_accumulation.cuh"
#include "cuda_nonbonded_force.cuh"
#include "geometry.h"

namespace {

__device__ __forceinline__ void emit_exact_packed_segment(
    int segment_begin,
    int segment_atoms,
    int segment_type_slot,

    int target_start,
    int target_size,

    int& entry_cursor,
    ExactEntry* __restrict__ exact_entries) {
    if (segment_atoms == 0) {
        return;
    }

    /*
     * X is the fixed, contiguous target group.
     * Y is the packed indirect atom list that rotates.
     */
    for (int x_offset = 0; x_offset < target_size; x_offset += 32) {
        const int x_len = min(32, target_size - x_offset);

        for (int y_offset = 0; y_offset < segment_atoms; y_offset += 32) {
            const int y_len = min(32, segment_atoms - y_offset);

            ExactEntry entry{};

            entry.x_start = target_start + x_offset;
            entry.x_len = x_len;

            entry.y_start = segment_begin + y_offset;
            entry.y_len = y_len;

            entry.y_type_slot = segment_type_slot;

            entry.diagonal = 0;
            entry.y_indirect = 1;

            exact_entries[entry_cursor++] = entry;
        }
    }
}

__device__ bool same_exact_energy_class(int slot1, int slot2, const uint8_t* category, const int* q_state) {
    return category[slot1] == category[slot2] && q_state[slot1] == q_state[slot2];
}

__device__ __forceinline__ void write_unique_lrf_component(LrfCoefficients& output, int phi, double value) {
    switch (phi) {
        case 0:
            output.phi0 = value;
            break;

        case 1:
            output.phi1[0] = value;
            break;

        case 2:
            output.phi1[1] = value;
            break;

        case 3:
            output.phi1[2] = value;
            break;

        case 4:
            output.phi2[0] = value;
            break;

        case 5:
            output.phi2[4] = value;
            break;

        case 6:
            output.phi2[8] = value;
            break;

        case 7:
            output.phi2[1] = value;
            output.phi2[3] = value;
            break;

        case 8:
            output.phi2[2] = value;
            output.phi2[6] = value;
            break;

        case 9:
            output.phi2[5] = value;
            output.phi2[7] = value;
            break;

        case 10:
            output.phi3[0] = value;
            break;

        case 11:
            output.phi3[13] = value;
            break;

        case 12:
            output.phi3[26] = value;
            break;

        /*
         * xxy:
         * (x,x,y), (x,y,x), (y,x,x)
         */
        case 13:
            output.phi3[1] = value;
            output.phi3[3] = value;
            output.phi3[9] = value;
            break;

        /*
         * xxz:
         * (x,x,z), (x,z,x), (z,x,x)
         */
        case 14:
            output.phi3[2] = value;
            output.phi3[6] = value;
            output.phi3[18] = value;
            break;

        /*
         * xyy:
         * (x,y,y), (y,x,y), (y,y,x)
         */
        case 15:
            output.phi3[4] = value;
            output.phi3[10] = value;
            output.phi3[12] = value;
            break;

        /*
         * yyz:
         * (y,y,z), (y,z,y), (z,y,y)
         */
        case 16:
            output.phi3[14] = value;
            output.phi3[16] = value;
            output.phi3[22] = value;
            break;

        /*
         * xzz:
         * (x,z,z), (z,x,z), (z,z,x)
         */
        case 17:
            output.phi3[8] = value;
            output.phi3[20] = value;
            output.phi3[24] = value;
            break;

        /*
         * yzz:
         * (y,z,z), (z,y,z), (z,z,y)
         */
        case 18:
            output.phi3[17] = value;
            output.phi3[23] = value;
            output.phi3[25] = value;
            break;

        /*
         * xyz: all six permutations.
         */
        case 19:
            output.phi3[5] = value;
            output.phi3[7] = value;
            output.phi3[11] = value;
            output.phi3[15] = value;
            output.phi3[19] = value;
            output.phi3[21] = value;
            break;
    }
}

template <int ORDER, int COMPONENT_COUNT>
__global__ void build_lrf_coefficients_order_csr_kernel(
    int n_group_ranges,

    const int* __restrict__ group_indices,
    const int* __restrict__ group_start_idx,
    const uint8_t* __restrict__ category,
    const int* __restrict__ atom_offsets,

    const double* __restrict__ dx,
    const double* __restrict__ dy,
    const double* __restrict__ dz,

    const double* __restrict__ q_r1,
    const double* __restrict__ q_r3,
    const double* __restrict__ q_r5,
    const double* __restrict__ q_r7,

    LrfCoefficients* __restrict__ coefficients) {
    constexpr int WARPS_PER_BLOCK = 4;
    constexpr unsigned FULL_MASK = 0xffffffffu;

    const int target_range = blockIdx.x;
    const int thread = threadIdx.x;
    const int lane = thread & 31;
    const int warp = thread >> 5;

    if (target_range >= n_group_ranges) {
        return;
    }

    const int target_start = group_start_idx[target_range];
    const uint8_t target_category = category[target_start];

    constexpr uint8_t P = static_cast<uint8_t>(AtomCategory::P);
    constexpr uint8_t W = static_cast<uint8_t>(AtomCategory::W);

    if (target_category != P && target_category != W) {
        return;
    }

    const int begin = atom_offsets[target_range];
    const int end = atom_offsets[target_range + 1];

    double sums[COMPONENT_COUNT];

#pragma unroll
    for (int component = 0; component < COMPONENT_COUNT; ++component) {
        sums[component] = 0.0;
    }

    for (int entry = begin + thread; entry < end; entry += blockDim.x) {
        if constexpr (ORDER == 0) {
            sums[0] += q_r1[entry];
        }

        if constexpr (ORDER == 1) {
            const double x = dx[entry];
            const double y = dy[entry];
            const double z = dz[entry];
            const double qr3 = q_r3[entry];

            sums[0] -= x * qr3;
            sums[1] -= y * qr3;
            sums[2] -= z * qr3;
        }

        if constexpr (ORDER == 2) {
            const double x = dx[entry];
            const double y = dy[entry];
            const double z = dz[entry];

            const double qr3 = q_r3[entry];
            const double qr5 = q_r5[entry];
            const double three_qr5 = 3.0 * qr5;

            sums[0] += x * x * three_qr5 - qr3;
            sums[1] += y * y * three_qr5 - qr3;
            sums[2] += z * z * three_qr5 - qr3;

            sums[3] += x * y * three_qr5;
            sums[4] += x * z * three_qr5;
            sums[5] += y * z * three_qr5;
        }

        if constexpr (ORDER == 3) {
            const double x = dx[entry];
            const double y = dy[entry];
            const double z = dz[entry];

            const double qr5 = q_r5[entry];
            const double qr7 = q_r7[entry];

            const double three_qr5 = 3.0 * qr5;
            const double nine_qr5 = 9.0 * qr5;
            const double fifteen_qr7 = 15.0 * qr7;

            const double xx = x * x;
            const double yy = y * y;
            const double zz = z * z;

            sums[0] += x * nine_qr5 - x * xx * fifteen_qr7;

            sums[1] += y * nine_qr5 - y * yy * fifteen_qr7;

            sums[2] += z * nine_qr5 - z * zz * fifteen_qr7;

            sums[3] += y * three_qr5 - xx * y * fifteen_qr7;

            sums[4] += z * three_qr5 - xx * z * fifteen_qr7;

            sums[5] += x * three_qr5 - x * yy * fifteen_qr7;

            sums[6] += z * three_qr5 - yy * z * fifteen_qr7;

            sums[7] += x * three_qr5 - x * zz * fifteen_qr7;

            sums[8] += y * three_qr5 - y * zz * fifteen_qr7;

            sums[9] -= x * y * z * fifteen_qr7;
        }
    }

#pragma unroll
    for (int component = 0; component < COMPONENT_COUNT; ++component) {
#pragma unroll
        for (int offset = 16; offset > 0; offset >>= 1) {
            sums[component] += __shfl_down_sync(FULL_MASK, sums[component], offset);
        }
    }

    __shared__ double warp_sums[COMPONENT_COUNT][WARPS_PER_BLOCK];

    if (lane == 0) {
#pragma unroll
        for (int component = 0; component < COMPONENT_COUNT; ++component) {
            warp_sums[component][warp] = sums[component];
        }
    }

    __syncthreads();

    if (warp != 0) {
        return;
    }

#pragma unroll
    for (int component = 0; component < COMPONENT_COUNT; ++component) {
        double block_sum = lane < WARPS_PER_BLOCK ? warp_sums[component][lane] : 0.0;

#pragma unroll
        for (int offset = 16; offset > 0; offset >>= 1) {
            block_sum += __shfl_down_sync(FULL_MASK, block_sum, offset);
        }

        if (lane == 0) {
            constexpr int first_component = ORDER == 0 ? 0 : ORDER == 1 ? 1
                                                         : ORDER == 2   ? 4
                                                                        : 10;

            write_unique_lrf_component(coefficients[group_indices[target_range]], first_component + component, block_sum);
        }
    }
}

__device__ int get_pair_index(
    int n,
    int group1,
    int group2) {
    const int x = min(group1, group2);
    const int y = max(group1, group2);

    // Row x starts after:
    // n + (n-1) + ... + (n-x+1)
    return x * n - (x * (x + 1)) / 2 + y;
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
    real_t inv_dis2 = static_cast<real_t>(1.0) / dis2;
    real_t inv_dis = sqrt(inv_dis2);

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

__device__ void shuffle(int& atom, uint8_t& atom_type, int& atom_state, real_t& atom_charge, vdw_atom_param_t& atom_vdw, real_t& atom_lambda, real_t3& atom_force, real_t3& atom_coord) {
    constexpr unsigned FULL_MASK = 0xFFFFFFFF;
    int src = ((threadIdx.x & 31) + 1) & 31;
    atom = __shfl_sync(FULL_MASK, atom, src);
    int tmp = atom_type;
    atom_type = static_cast<uint8_t>(__shfl_sync(FULL_MASK, tmp, src));
    atom_state = __shfl_sync(FULL_MASK, atom_state, src);
    atom_charge = __shfl_sync(FULL_MASK, atom_charge, src);

    atom_vdw.aii_normal = __shfl_sync(FULL_MASK, atom_vdw.aii_normal, src);
    atom_vdw.bii_normal = __shfl_sync(FULL_MASK, atom_vdw.bii_normal, src);
    atom_vdw.aii_14 = __shfl_sync(FULL_MASK, atom_vdw.aii_14, src);
    atom_vdw.bii_14 = __shfl_sync(FULL_MASK, atom_vdw.bii_14, src);

    atom_lambda = __shfl_sync(FULL_MASK, atom_lambda, src);
    atom_force.x = __shfl_sync(FULL_MASK, atom_force.x, src);
    atom_force.y = __shfl_sync(FULL_MASK, atom_force.y, src);
    atom_force.z = __shfl_sync(FULL_MASK, atom_force.z, src);
    atom_coord.x = __shfl_sync(FULL_MASK, atom_coord.x, src);
    atom_coord.y = __shfl_sync(FULL_MASK, atom_coord.y, src);
    atom_coord.z = __shfl_sync(FULL_MASK, atom_coord.z, src);
}

__device__ void nonbonded_force_calculation(
    int x_idx,
    int y_idx,
    bool is_diag,
    int base_x,
    int base_y,

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
    energy_accum_t* e

) {
    int lane = threadIdx.x & 31;
    const int atom1 = x_idx == -1 ? -1 : atom_idx[x_idx];

    const auto& atom1_type = atom1 == -1 ? static_cast<uint8_t>(AtomCategory::INVALID) : category[x_idx];
    const int atom1_state = atom1 == -1 ? -1 : q_state[x_idx];
    const real_t atom1_charge = atom1 == -1 ? 0 : atom_charge[x_idx];
    const vdw_atom_param_t atom1_vdw = atom1 == -1 ? vdw_atom_param_t{0, 0, 0, 0} : atom_vdw[x_idx];
    const real_t atom1_lambda = atom1 == -1 ? 0 : atom_lambdas[x_idx];
    real_t3 atom1_coord = atom1 == -1 ? real_t3{0, 0, 0} : real_t3{cx[x_idx], cy[x_idx], cz[x_idx]};
    real_t3 atom1_force = {0, 0, 0};

    int atom2 = y_idx == -1 ? -1 : atom_idx[y_idx];
    uint8_t atom2_type = atom2 == -1 ? static_cast<uint8_t>(AtomCategory::INVALID) : category[y_idx];
    int atom2_state = atom2 == -1 ? -1 : q_state[y_idx];
    real_t atom2_charge = atom2 == -1 ? 0 : atom_charge[y_idx];
    vdw_atom_param_t atom2_vdw = atom2 == -1 ? vdw_atom_param_t{0, 0, 0, 0} : atom_vdw[y_idx];
    real_t atom2_lambda = atom2 == -1 ? 0 : atom_lambdas[y_idx];
    real_t3 atom2_coord = atom2 == -1 ? real_t3{0, 0, 0} : real_t3{cx[y_idx], cy[y_idx], cz[y_idx]};
    real_t3 atom2_force = {0, 0, 0};

    real_t local_e_coul = 0, local_e_vdw = 0;
    for (int i = 0; i < 32; i++) {
        if (!is_diag || atom1 < atom2) {
            compute_pair(atom1, atom1_type, atom1_state, atom1_charge, atom1_vdw, atom1_lambda, atom1_coord,
                         atom2, atom2_type, atom2_state, atom2_charge, atom2_vdw, atom2_lambda, atom2_coord,
                         n_atoms_solute, LJ_matrix,
                         el14_scale, coulomb_constant, vdw_rule,
                         atom1_force, atom2_force,
                         local_e_coul, local_e_vdw);
        }
        shuffle(atom2, atom2_type, atom2_state, atom2_charge, atom2_vdw, atom2_lambda, atom2_force, atom2_coord);
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
    const unsigned mask = 0xffffffffu;
    for (int offset = 16; offset > 0; offset >>= 1) {
        local_e_coul += __shfl_down_sync(mask, local_e_coul, offset);
        local_e_vdw += __shfl_down_sync(mask, local_e_vdw, offset);
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

__global__ void count_exact_atom_tiles_kernel(
    int n_group_ranges,

    const uint8_t* group_pair_modes,
    const int* group_start_idx,
    const int* group_sizes,
    const uint8_t* category,
    const int* q_state,

    int* atom_degrees,
    int* entry_degrees) {
    const int target_range = blockIdx.x * blockDim.x + threadIdx.x;

    if (target_range >= n_group_ranges) {
        return;
    }

    const int target_size = group_sizes[target_range];

    /*
     * The target group is the fixed, contiguous X axis.
     */
    const int target_x_tiles = (target_size + 31) / 32;

    int atom_count = 0;
    int entry_count = 0;

    int segment_atoms = 0;
    int segment_type_slot = -1;

    /*
     * source_range < target_range:
     * only cross-group atoms are packed into the indirect Y list.
     */
    for (int source_range = 0; source_range < target_range; ++source_range) {
        const int pair_index = get_pair_index(n_group_ranges, source_range, target_range);

        if (group_pair_modes[pair_index] != GROUP_PAIR_EXACT) {
            continue;
        }

        const int source_start = group_start_idx[source_range];

        const int source_size = group_sizes[source_range];

        /*
         * Finish the current packed Y segment when its
         * category or Q state changes.
         */
        if (segment_atoms != 0 && !same_exact_energy_class(segment_type_slot, source_start, category, q_state)) {
            const int segment_y_tiles = (segment_atoms + 31) / 32;

            entry_count += target_x_tiles * segment_y_tiles;

            segment_atoms = 0;
            segment_type_slot = -1;
        }

        if (segment_atoms == 0) {
            segment_type_slot = source_start;
        }

        segment_atoms += source_size;
        atom_count += source_size;
    }

    /*
     * Count the final packed Y segment.
     */
    if (segment_atoms != 0) {
        const int segment_y_tiles = (segment_atoms + 31) / 32;

        entry_count += target_x_tiles * segment_y_tiles;
    }

    /*
     * Count the target group's separate diagonal entries.
     * These entries do not use the packed Y list.
     */
    const int diagonal_pair = get_pair_index(n_group_ranges, target_range, target_range);

    if (group_pair_modes[diagonal_pair] == GROUP_PAIR_EXACT) {
        entry_count += target_x_tiles * (target_x_tiles + 1) / 2;
    }

    atom_degrees[target_range] = atom_count;
    entry_degrees[target_range] = entry_count;
}

__global__ void count_lrf_atom_degree_kernel(int n_group_ranges, const uint8_t* group_pair_modes, const int* group_sizes, int* degrees) {
    const int target_range = blockIdx.x * blockDim.x + threadIdx.x;

    if (target_range >= n_group_ranges) {
        return;
    }

    int atom_count = 0;
    for (int source_range = 0; source_range < n_group_ranges; source_range++) {
        const int pair_index = get_pair_index(n_group_ranges, target_range, source_range);
        if (group_pair_modes[pair_index] == GROUP_PAIR_LRF) {
            atom_count += group_sizes[source_range];
        }
    }
    degrees[target_range] = atom_count;
}

__global__ void fill_lrf_atom_csr_kernel(int n_group_ranges, const uint8_t* group_pair_modes, const int* group_start_idx, const int* group_sizes, const int* offsets, int* source_atom_slots) {
    const int target_range = blockIdx.x * blockDim.x + threadIdx.x;

    if (target_range >= n_group_ranges) {
        return;
    }

    int output = offsets[target_range];

    for (int source_range = 0; source_range < n_group_ranges; source_range++) {
        const int pair_index = get_pair_index(n_group_ranges, target_range, source_range);
        if (group_pair_modes[pair_index] != GROUP_PAIR_LRF) {
            continue;
        }
        const int source_start = group_start_idx[source_range];
        const int source_size = group_sizes[source_range];
        for (int local_atom = 0; local_atom < source_size; local_atom++) {
            source_atom_slots[output++] = source_start + local_atom;
        }
    }
}

__global__ void calc_lrf_kernel(
    int n_slots,

    const int* atom_idx,
    const int* atom_to_group,
    const uint8_t* category,
    const real_t* atom_charge,

    const real_t* cx,
    const real_t* cy,
    const real_t* cz,

    const LrfCoefficients* coefficients,
    double coulomb_constant,

    dvel_t* dvelocities,
    energy_accum_t* energy) {
    const int slot = blockIdx.x * blockDim.x + threadIdx.x;

    if (slot >= n_slots) {
        return;
    }

    constexpr uint8_t P = static_cast<uint8_t>(AtomCategory::P);

    constexpr uint8_t W = static_cast<uint8_t>(AtomCategory::W);

    const uint8_t atom_category = category[slot];

    if (atom_category != P && atom_category != W) {
        return;
    }

    const int atom = atom_idx[slot];

    if (atom < 0) {
        return;
    }

    const int group = atom_to_group[atom];

    if (group < 0) {
        return;
    }

    const LrfCoefficients& lrf = coefficients[group];
    const double d[3] = {
        lrf.center.x - static_cast<double>(cx[slot]),
        lrf.center.y - static_cast<double>(cy[slot]),
        lrf.center.z - static_cast<double>(cz[slot])};

    double potential = lrf.phi0;

    for (int a = 0; a < 3; ++a) {
        potential += lrf.phi1[a] * d[a];
    }

    for (int a = 0; a < 3; ++a) {
        for (int b = 0; b < 3; ++b) {
            potential += 0.5 * lrf.phi2[a * 3 + b] * d[a] * d[b];
        }
    }

    double df[3] = {lrf.phi1[0], lrf.phi1[1], lrf.phi1[2]};

    for (int a = 0; a < 3; ++a) {
        for (int b = 0; b < 3; ++b) {
            df[a] += lrf.phi2[a * 3 + b] * d[b];
        }
    }

    for (int a = 0; a < 3; ++a) {
        for (int b = 0; b < 3; ++b) {
            for (int c = 0; c < 3; ++c) {
                const int index = (a * 3 + b) * 3 + c;
                df[a] += 0.5 * lrf.phi3[index] * d[b] * d[c];
            }
        }
    }

    const double charge = static_cast<double>(atom_charge[slot]);

    const double energy_value = 0.5 * coulomb_constant * charge * potential;

    const double force_scale = -coulomb_constant * charge;

    atomic_add_force(&dvelocities[atom].x, force_scale * df[0]);

    atomic_add_force(&dvelocities[atom].y, force_scale * df[1]);

    atomic_add_force(&dvelocities[atom].z, force_scale * df[2]);

    atomic_add_energy(&energy[E_LRF], energy_value);
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

    auto [tile_x, tile_y] = get_tile_idx(block_num, tile);

    const int base_x = tile_x << 5;
    const int base_y = tile_y << 5;

    int x_idx = base_x + lane;
    int y_idx = base_y + lane;
    x_idx = x_idx < sz ? x_idx : -1;
    y_idx = y_idx < sz ? y_idx : -1;

    nonbonded_force_calculation(x_idx, y_idx, tile_x == tile_y, base_x, base_y, n_states, n_atoms_solute,
                                atom_idx, category, q_state, atom_lambdas, atom_charge, atom_vdw, LJ_matrix, el14_scale, coulomb_constant, vdw_rule, cx, cy, cz, dvelocities, e);
}

__global__ void classify_group_pairs_by_switch_kernel(
    int n_groups_ranges,

    double solute_solute_cutoff2,
    double solute_solvent_cutoff2,
    double solvent_solvent_cutoff2,
    double rcq2,
    double lrf_cutoff2,

    const coord_t solute_center,
    const int* group_start_idx,
    const int* atom_idx,
    const uint8_t* category,
    const int* q_state,
    const coord_t* coords,

    uint8_t* group_pair_modes) {
    const int pair_index = blockIdx.x * blockDim.x + threadIdx.x;
    const int total_pairs = n_groups_ranges * (n_groups_ranges + 1) / 2;

    if (pair_index >= total_pairs) {
        return;
    }

    const int2 pair = get_tile_idx(n_groups_ranges, pair_index);
    const int group1 = pair.x;
    const int group2 = pair.y;

    group_pair_modes[pair_index] = GROUP_PAIR_IGNORE;
    const int switch_atom1 = group_start_idx[group1];
    const int switch_atom2 = group_start_idx[group2];

    const uint8_t category1 = category[switch_atom1];
    const uint8_t category2 = category[switch_atom2];

    constexpr uint8_t P = static_cast<uint8_t>(AtomCategory::P);
    constexpr uint8_t Q = static_cast<uint8_t>(AtomCategory::Q);
    constexpr uint8_t W = static_cast<uint8_t>(AtomCategory::W);

    const bool group1_is_q = category1 == Q;
    const bool group2_is_q = category2 == Q;

    if (group1_is_q && group2_is_q) {
        const int state1 = q_state[switch_atom1];
        const int state2 = q_state[switch_atom2];
        if (state1 == state2) {
            group_pair_modes[pair_index] = GROUP_PAIR_EXACT;
        }
        return;
    }

    if (group1_is_q || group2_is_q) {
        // Q-P or Q-W
        const int environment_switch_atom = group1_is_q ? switch_atom2 : switch_atom1;
        const double environment_distance2 = norm2(coords[atom_idx[environment_switch_atom]] - solute_center);

        if (environment_distance2 <= rcq2) {
            group_pair_modes[pair_index] = GROUP_PAIR_EXACT;
        }
        return;
    }

    // P-P, P-W, or W-W
    const double group_distance2 = norm2(coords[atom_idx[switch_atom1]] - coords[atom_idx[switch_atom2]]);
    double normal_cutoff2;
    if (category1 == P && category2 == P) {
        normal_cutoff2 = solute_solute_cutoff2;
    } else if (category1 == W && category2 == W) {
        normal_cutoff2 = solvent_solvent_cutoff2;
    } else {
        normal_cutoff2 = solute_solvent_cutoff2;
    }

    if (group_distance2 <= normal_cutoff2) {
        group_pair_modes[pair_index] = GROUP_PAIR_EXACT;
    } else if (group_distance2 <= lrf_cutoff2) {
        group_pair_modes[pair_index] = GROUP_PAIR_LRF;
    }
}

__global__ void exact_tiles_nonbonded_force_kernel(
    int n_exact_tiles,
    const ExactEntry* exact_entries,
    const int* exact_source_atom_slots,

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
    const int lane = threadIdx.x & 31;
    const int warp_in_block = threadIdx.x >> 5;
    const int warps_per_block = blockDim.x >> 5;

    const int tile_index = blockIdx.x * warps_per_block + warp_in_block;

    if (tile_index >= n_exact_tiles) {
        return;
    }

    const ExactEntry tile = exact_entries[tile_index];

    /*
     * Each lane loads one X atom and one Y atom.
     * Invalid lanes use atom == -1 and still participate in all shuffles.
     */
    const int x_idx = lane < tile.x_len ? tile.x_start + lane : -1;
    const int y_idx = lane < tile.y_len ? tile.y_indirect == 1 ? exact_source_atom_slots[tile.y_start + lane] : tile.y_start + lane : -1;

    nonbonded_force_calculation(x_idx, y_idx, tile.diagonal == 1, tile.x_start, tile.y_type_slot, n_states, n_atoms_solute, atom_idx,
                                category, q_state, atom_lambdas, atom_charge, atom_vdw, LJ_matrix, el14_scale, coulomb_constant, vdw_rule, cx, cy, cz, dvelocities, e);
}

__global__ void compute_lrf_centers_kernel(
    int n_group_ranges,

    const int* group_indices,
    const int* group_start_idx,
    const int* group_sizes,

    const uint8_t* category,

    const real_t* cx,
    const real_t* cy,
    const real_t* cz,

    LrfCoefficients* coefficients) {
    const int lane = threadIdx.x & 31;
    const int warp_in_block = threadIdx.x >> 5;
    const int warps_per_block = blockDim.x >> 5;

    const int range = blockIdx.x * warps_per_block + warp_in_block;

    if (range >= n_group_ranges) {
        return;
    }

    const int start = group_start_idx[range];
    const int size = group_sizes[range];

    constexpr uint8_t P = static_cast<uint8_t>(AtomCategory::P);
    constexpr uint8_t W = static_cast<uint8_t>(AtomCategory::W);

    const uint8_t type = category[start];

    if (type != P && type != W) {
        return;
    }

    double sum_x = 0;
    double sum_y = 0;
    double sum_z = 0;

    for (int i = lane; i < size; i += 32) {
        const int slot = start + i;

        sum_x += static_cast<double>(cx[slot]);
        sum_y += static_cast<double>(cy[slot]);
        sum_z += static_cast<double>(cz[slot]);
    }

    constexpr unsigned MASK = 0xffffffffu;

    for (int offset = 16; offset > 0; offset >>= 1) {
        sum_x += __shfl_down_sync(MASK, sum_x, offset);
        sum_y += __shfl_down_sync(MASK, sum_y, offset);
        sum_z += __shfl_down_sync(MASK, sum_z, offset);
    }

    if (lane == 0) {
        const int original_group = group_indices[range];

        const double inv_size = 1.0 / static_cast<double>(size);

        coefficients[original_group].center = {
            sum_x * inv_size,
            sum_y * inv_size,
            sum_z * inv_size};
    }
}

__global__ void fill_exact_atom_tiles_kernel(
    int n_group_ranges,

    const uint8_t* group_pair_modes,
    const int* group_start_idx,
    const int* group_sizes,
    const uint8_t* category,
    const int* q_state,

    const int* atom_offsets,
    const int* entry_offsets,

    int* source_atom_slots,
    ExactEntry* exact_entries) {
    const int target_range = blockIdx.x * blockDim.x + threadIdx.x;
    if (target_range >= n_group_ranges) {
        return;
    }

    const int target_start = group_start_idx[target_range];
    const int target_size = group_sizes[target_range];

    int atom_cursor = atom_offsets[target_range];
    int entry_cursor = entry_offsets[target_range];

    int segment_begin = atom_cursor;
    int segment_atoms = 0;
    int segment_type_slot = -1;

    for (int source_range = 0; source_range < target_range; source_range++) {
        const int pair_index = get_pair_index(n_group_ranges, source_range, target_range);

        if (group_pair_modes[pair_index] != GROUP_PAIR_EXACT) {
            continue;
        }

        const int source_start = group_start_idx[source_range];

        const int source_size = group_sizes[source_range];

        if (segment_atoms != 0 && !same_exact_energy_class(segment_type_slot, source_start, category, q_state)) {
            emit_exact_packed_segment(segment_begin, segment_atoms, segment_type_slot, target_start, target_size, entry_cursor, exact_entries);
            segment_begin = atom_cursor;
            segment_atoms = 0;
            segment_type_slot = -1;
        }

        if (segment_atoms == 0) {
            segment_begin = atom_cursor;
            segment_type_slot = source_start;
        }

        for (int local_atom = 0; local_atom < source_size; local_atom++) {
            source_atom_slots[atom_cursor++] = source_start + local_atom;
            segment_atoms++;
        }
    }

    emit_exact_packed_segment(segment_begin, segment_atoms, segment_type_slot, target_start, target_size, entry_cursor, exact_entries);

    const int diagonal_pair = get_pair_index(n_group_ranges, target_range, target_range);

    if (group_pair_modes[diagonal_pair] == GROUP_PAIR_EXACT) {
        const int n = (target_size + 31) / 32;

        for (int ix = 0; ix < n; ix++) {
            const int x_offset = ix * 32;
            const int x_len = min(32, target_size - x_offset);

            for (int iy = ix; iy < n; ++iy) {
                const int y_offset = iy * 32;
                const int y_len = min(32, target_size - y_offset);

                ExactEntry entry{};

                entry.x_start = target_start + x_offset;

                entry.x_len = x_len;
                entry.y_start = target_start + y_offset;
                entry.y_type_slot = entry.y_start;

                entry.y_len = y_len;

                entry.diagonal = (ix == iy);
                entry.y_indirect = 0;

                exact_entries[entry_cursor++] = entry;
            }
        }
    }
}

constexpr int LRF_UNIQUE_COMPONENTS = 20;
constexpr int LRF_COEFFICIENT_THREADS = 256;

__global__ void build_lrf_coefficients_dense_kernel(
    int n_group_ranges,
    int n_slots,

    const int* __restrict__ group_indices,
    const int* __restrict__ group_start_idx,
    const uint8_t* __restrict__ category,
    const uint8_t* __restrict__ group_pair_modes,
    const int* __restrict__ slot_to_group_range,

    const real_t* __restrict__ atom_charge,
    const real_t* __restrict__ cx,
    const real_t* __restrict__ cy,
    const real_t* __restrict__ cz,

    LrfCoefficients* __restrict__ coefficients) {
    constexpr int WARPS_PER_BLOCK = LRF_COEFFICIENT_THREADS / 32;

    constexpr unsigned FULL_MASK = 0xffffffffu;

    /*
     * Exactly one block processes one target range.
     */
    const int target_range = blockIdx.x;

    if (target_range >= n_group_ranges) {
        return;
    }

    const int target_start = group_start_idx[target_range];

    const uint8_t target_category = category[target_start];

    constexpr uint8_t P = static_cast<uint8_t>(AtomCategory::P);

    constexpr uint8_t W = static_cast<uint8_t>(AtomCategory::W);

    /*
     * LRF is only evaluated for P and W atoms.
     * This condition is uniform across the whole block.
     */
    if (target_category != P && target_category != W) {
        return;
    }

    const int target_group = group_indices[target_range];

    const coord_t center = coefficients[target_group].center;

    const int thread = threadIdx.x;
    const int lane = thread & 31;
    const int warp = thread >> 5;

    /*
     * The 20 independent components are:
     *
     *  0: phi0
     *
     *  1: x
     *  2: y
     *  3: z
     *
     *  4: xx
     *  5: yy
     *  6: zz
     *  7: xy
     *  8: xz
     *  9: yz
     *
     * 10: xxx
     * 11: yyy
     * 12: zzz
     * 13: xxy
     * 14: xxz
     * 15: xyy
     * 16: yyz
     * 17: xzz
     * 18: yzz
     * 19: xyz
     */
    double sums[LRF_UNIQUE_COMPONENTS];

#pragma unroll
    for (int component = 0; component < LRF_UNIQUE_COMPONENTS; ++component) {
        sums[component] = 0.0;
    }

    /*
     * Threads collectively scan the complete contiguous source-slot
     * array. Padding slots have source_range == -1.
     */
    for (int source_slot = thread; source_slot < n_slots; source_slot += blockDim.x) {
        const int source_range = slot_to_group_range[source_slot];

        if (source_range < 0) {
            continue;
        }

        const int pair_index = get_pair_index(n_group_ranges, target_range, source_range);

        /*
         * This preserves the original CSR semantics exactly:
         *
         * - GROUP_PAIR_LRF is included.
         * - GROUP_PAIR_EXACT is excluded.
         * - GROUP_PAIR_IGNORE is excluded.
         */
        if (group_pair_modes[pair_index] != GROUP_PAIR_LRF) {
            continue;
        }

        const double charge = static_cast<double>(atom_charge[source_slot]);

        const double x = static_cast<double>(cx[source_slot]) - center.x;

        const double y = static_cast<double>(cy[source_slot]) - center.y;

        const double z = static_cast<double>(cz[source_slot]) - center.z;

        const double r2 = x * x + y * y + z * z;

        /*
         * A self relationship should not be LRF, but guard against
         * singular values in case classification changes later.
         */
        if (r2 == 0.0) {
            continue;
        }

        const double inv_r = rsqrt(r2);
        const double inv_r2 = 1.0 / r2;

        const double inv_r3 = inv_r * inv_r2;

        const double inv_r5 = inv_r3 * inv_r2;

        const double inv_r7 = inv_r5 * inv_r2;

        const double q_r1 = charge * inv_r;

        const double q_r3 = charge * inv_r3;

        const double q_r5 = charge * inv_r5;

        const double q_r7 = charge * inv_r7;

        const double three_q_r5 = 3.0 * q_r5;

        const double nine_q_r5 = 9.0 * q_r5;

        const double fifteen_q_r7 = 15.0 * q_r7;

        const double xx = x * x;
        const double yy = y * y;
        const double zz = z * z;

        /*
         * phi0
         */
        sums[0] += q_r1;

        /*
         * phi1
         */
        sums[1] -= x * q_r3;
        sums[2] -= y * q_r3;
        sums[3] -= z * q_r3;

        /*
         * Unique phi2 components.
         */
        sums[4] += xx * three_q_r5 - q_r3;

        sums[5] += yy * three_q_r5 - q_r3;

        sums[6] += zz * three_q_r5 - q_r3;

        sums[7] += x * y * three_q_r5;

        sums[8] += x * z * three_q_r5;

        sums[9] += y * z * three_q_r5;

        /*
         * Unique phi3 components.
         */
        sums[10] += x * nine_q_r5 - x * xx * fifteen_q_r7;

        sums[11] += y * nine_q_r5 - y * yy * fifteen_q_r7;

        sums[12] += z * nine_q_r5 - z * zz * fifteen_q_r7;

        sums[13] += y * three_q_r5 - xx * y * fifteen_q_r7;

        sums[14] += z * three_q_r5 - xx * z * fifteen_q_r7;

        sums[15] += x * three_q_r5 - x * yy * fifteen_q_r7;

        sums[16] += z * three_q_r5 - yy * z * fifteen_q_r7;

        sums[17] += x * three_q_r5 - x * zz * fifteen_q_r7;

        sums[18] += y * three_q_r5 - y * zz * fifteen_q_r7;

        sums[19] -= x * y * z * fifteen_q_r7;
    }

    /*
     * Reduce each component within every warp.
     *
     * This is performed once after all source atoms have been
     * processed, instead of once per source chunk.
     */
#pragma unroll
    for (int component = 0; component < LRF_UNIQUE_COMPONENTS; ++component) {
#pragma unroll
        for (int offset = 16; offset > 0; offset >>= 1) {
            sums[component] += __shfl_down_sync(FULL_MASK, sums[component], offset);
        }
    }

    __shared__ double warp_sums[LRF_UNIQUE_COMPONENTS][WARPS_PER_BLOCK];

    if (lane == 0) {
#pragma unroll
        for (int component = 0; component < LRF_UNIQUE_COMPONENTS; ++component) {
            warp_sums[component][warp] = sums[component];
        }
    }

    __syncthreads();

    /*
     * Warp zero reduces the eight warp results.
     */
    if (warp == 0) {
#pragma unroll
        for (int component = 0; component < LRF_UNIQUE_COMPONENTS; ++component) {
            double value = lane < WARPS_PER_BLOCK ? warp_sums[component][lane] : 0.0;

#pragma unroll
            for (int offset = 16; offset > 0; offset >>= 1) {
                value += __shfl_down_sync(FULL_MASK, value, offset);
            }

            if (lane == 0) {
                write_unique_lrf_component(coefficients[target_group], component, value);
            }
        }
    }
}

}  // namespace

void CudaNonbondedForce::init_backend(Context& ctx) {
    // Buffers are indexed by combined-list position [0, n_total), which exceeds
    // n_atoms because Q atoms are duplicated per FEP state and the list is padded.
    coord_x_ = std::make_unique<HostDeviceBuffer<real_t>>(data_.n_total);
    coord_y_ = std::make_unique<HostDeviceBuffer<real_t>>(data_.n_total);
    coord_z_ = std::make_unique<HostDeviceBuffer<real_t>>(data_.n_total);

    const int n_group_ranges = data_.group_indices->length;

    const int max_group_pairs = n_group_ranges * (n_group_ranges + 1) / 2;
    int max_exact_tiles = 0;

    const int* group_sizes = data_.group_sizes->cpu_data_p;

    for (int group1 = 0; group1 < n_group_ranges; group1++) {
        int nx = (group_sizes[group1] + 31) / 32;
        for (int group2 = group1; group2 < n_group_ranges; group2++) {
            int ny = (group_sizes[group2] + 31) / 32;
            if (group1 == group2) {
                max_exact_tiles += nx * (nx + 1) / 2;
            } else {
                max_exact_tiles += nx * ny;
            }
        }
    }

    exact_tile_capacity_ = max_exact_tiles;
    lrf_pair_capacity_ = max_group_pairs;

    group_pair_modes_ = std::make_unique<HostDeviceBuffer<uint8_t>>(max_group_pairs, false, true);

    exact_tiles_ = std::make_unique<HostDeviceBuffer<ExactEntry>>(exact_tile_capacity_, false, true);

    lrf_group_pairs_ = std::make_unique<HostDeviceBuffer<LrfPairEntry>>(lrf_pair_capacity_, false, true);

    exact_tile_count_ = std::make_unique<HostDeviceBuffer<int>>(1, true, true);

    lrf_pair_count_ = std::make_unique<HostDeviceBuffer<int>>(1, true, true);

    list_overflow_ = std::make_unique<HostDeviceBuffer<int>>(1, true, true);

    lrf_coefficients_ = std::make_unique<HostDeviceBuffer<LrfCoefficients>>(ctx.charge_group_config.charge_groups.size(), false, true);

    lrf_atom_degrees_ = std::make_unique<HostDeviceBuffer<int>>(n_group_ranges + 1, false, true);

    lrf_atom_offsets_ = std::make_unique<HostDeviceBuffer<int>>(n_group_ranges + 1, true, true);

    check_cuda(cub::DeviceScan::ExclusiveSum(
        nullptr,
        lrf_scan_temp_bytes_,
        lrf_atom_degrees_->gpu_data_p,
        lrf_atom_offsets_->gpu_data_p,
        n_group_ranges + 1));

    lrf_scan_temp_ = std::make_unique<HostDeviceBuffer<unsigned char>>(lrf_scan_temp_bytes_, false, true);

    exact_atom_degrees_ = std::make_unique<HostDeviceBuffer<int>>(n_group_ranges + 1, false, true);

    exact_atom_offsets_ = std::make_unique<HostDeviceBuffer<int>>(n_group_ranges + 1, true, true);

    exact_entry_degrees_ = std::make_unique<HostDeviceBuffer<int>>(n_group_ranges + 1, false, true);

    exact_entry_offsets_ = std::make_unique<HostDeviceBuffer<int>>(n_group_ranges + 1, true, true);

    std::vector<int> slot_to_group_range(data_.n_total, -1);

    for (int range = 0; range < n_group_ranges; ++range) {
        const int begin = data_.group_start_idx->cpu_data_p[range];

        const int size = data_.group_sizes->cpu_data_p[range];

        for (int local_atom = 0; local_atom < size; ++local_atom) {
            slot_to_group_range[begin + local_atom] = range;
        }
    }

    lrf_slot_to_group_range_ = HostDeviceBuffer<int>::from_vector(slot_to_group_range, ctx.command_info.requested_gpu);
}

void CudaNonbondedForce::calc_all_direct_pairs(Context& ctx) {
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
                                           coord_x_->gpu_data_p, coord_y_->gpu_data_p, coord_z_->gpu_data_p, ctx.dvelocities->gpu_data_p, ctx.energy.device());
}

void CudaNonbondedForce::init_calculation_groups_by_switch(Context& ctx) {
    exact_tile_count_->zero();
    lrf_pair_count_->zero();
    list_overflow_->zero();

    const double solute_solute_cutoff2 = ctx.md.solute_solute * ctx.md.solute_solute;
    const double solute_solvent_cutoff2 = ctx.md.solute_solvent * ctx.md.solute_solvent;
    const double solvent_solvent_cutoff2 = ctx.md.solvent_solvent * ctx.md.solvent_solvent;
    const double rcq2 = ctx.md.q_atom * ctx.md.q_atom;
    const double lrf_cutoff2 = ctx.md.lrf_cutoff * ctx.md.lrf_cutoff;

    const int thread_num = 256;
    const int n_group_ranges = data_.group_indices->length;
    int total_pairs = n_group_ranges * (n_group_ranges + 1) >> 1;
    int grid_sz = (total_pairs + thread_num - 1) / thread_num;

    dim3 grid = dim3(grid_sz);
    classify_group_pairs_by_switch_kernel<<<grid, thread_num>>>(n_group_ranges,
                                                                solute_solute_cutoff2,
                                                                solute_solvent_cutoff2,
                                                                solvent_solvent_cutoff2,
                                                                rcq2,
                                                                lrf_cutoff2,
                                                                ctx.topo.solute_center,
                                                                data_.group_start_idx->gpu_data_p,
                                                                data_.atom_idx->gpu_data_p,
                                                                data_.category->gpu_data_p,
                                                                data_.q_state->gpu_data_p,
                                                                ctx.coords->gpu_data_p,
                                                                group_pair_modes_->gpu_data_p);
    check_cuda(cudaGetLastError());

    build_exact_atom_tiles(ctx);
}

void CudaNonbondedForce::init_calculation_groups_by_all_atoms(Context& ctx) {
    throw std::runtime_error("CUDA LRF with iuse_switch_atom == 0 is not implemented yet");
}

void CudaNonbondedForce::init_calculation_groups(Context& ctx) {
    const auto& config = ctx.charge_group_config;
    if (config.iuse_switch_atom == 1) {
        // Use groups.iswitch to check the distance
        init_calculation_groups_by_switch(ctx);

    } else {
        // Should use every atoms to check the distance
        // init_calculation_groups_by_all_atoms(ctx);
        // todo: now alwasys use switch to test
        init_calculation_groups_by_switch(ctx);
    }
    // build_lrf_atom_csr(ctx);
}

void CudaNonbondedForce::calc_exact_tiles(Context& ctx) {
    if (n_exact_tiles_ <= 0) {
        return;
    }

    constexpr int thread_num = 256;
    constexpr int warps_per_block = thread_num / 32;

    const int grid_sz = (n_exact_tiles_ + warps_per_block - 1) / warps_per_block;

    exact_tiles_nonbonded_force_kernel<<<grid_sz, thread_num>>>(
        n_exact_tiles_,
        exact_tiles_->gpu_data_p,

        exact_source_atom_slots_->gpu_data_p,

        ctx.n_lambdas(),
        ctx.n_atoms_solute,

        data_.atom_idx->gpu_data_p,
        data_.category->gpu_data_p,
        data_.q_state->gpu_data_p,
        data_.atom_lambdas->gpu_data_p,
        data_.atom_charge->gpu_data_p,
        data_.atom_vdw->gpu_data_p,

        ctx.LJ_matrix->gpu_data_p,

        static_cast<real_t>(ctx.topo.el14_scale),

        static_cast<real_t>(ctx.topo.coulomb_constant),

        ctx.topo.vdw_rule,

        coord_x_->gpu_data_p,
        coord_y_->gpu_data_p,
        coord_z_->gpu_data_p,

        ctx.dvelocities->gpu_data_p,
        ctx.energy.device());

    check_cuda(cudaGetLastError());
}

void CudaNonbondedForce::init_lrf_coefficients(Context& ctx) {
    lrf_coefficients_->zero();

    const int n_group_ranges = static_cast<int>(data_.group_indices->length);

    if (n_group_ranges <= 0) {
        return;
    }

    /*
     * Calculate target-group centers first.
     */
    constexpr int center_threads = 256;
    constexpr int center_warps_per_block = center_threads / 32;

    const int center_grid = (n_group_ranges + center_warps_per_block - 1) / center_warps_per_block;

    compute_lrf_centers_kernel<<<center_grid, center_threads>>>(
        n_group_ranges,

        data_.group_indices->gpu_data_p,
        data_.group_start_idx->gpu_data_p,
        data_.group_sizes->gpu_data_p,
        data_.category->gpu_data_p,

        coord_x_->gpu_data_p,
        coord_y_->gpu_data_p,
        coord_z_->gpu_data_p,

        lrf_coefficients_->gpu_data_p);

    check_cuda(cudaGetLastError());

    /*
     * Exactly one block per target range.
     */
    build_lrf_coefficients_dense_kernel<<<n_group_ranges, LRF_COEFFICIENT_THREADS>>>(
        n_group_ranges,
        data_.n_total,

        data_.group_indices->gpu_data_p,
        data_.group_start_idx->gpu_data_p,
        data_.category->gpu_data_p,

        group_pair_modes_->gpu_data_p,

        lrf_slot_to_group_range_
            ->gpu_data_p,

        data_.atom_charge->gpu_data_p,

        coord_x_->gpu_data_p,
        coord_y_->gpu_data_p,
        coord_z_->gpu_data_p,

        lrf_coefficients_->gpu_data_p);

    check_cuda(cudaGetLastError());
}

void CudaNonbondedForce::calc_lrf(Context& ctx) {
    const int n_slots = data_.n_total;

    if (n_slots <= 0) {
        return;
    }

    constexpr int thread_num = 256;

    const int grid_sz = (n_slots + thread_num - 1) / thread_num;

    calc_lrf_kernel<<<grid_sz, thread_num>>>(
        n_slots,

        data_.atom_idx->gpu_data_p,
        data_.atom_to_group->gpu_data_p,
        data_.category->gpu_data_p,
        data_.atom_charge->gpu_data_p,

        coord_x_->gpu_data_p,
        coord_y_->gpu_data_p,
        coord_z_->gpu_data_p,

        lrf_coefficients_->gpu_data_p,

        ctx.topo.coulomb_constant,

        ctx.dvelocities->gpu_data_p,
        ctx.energy.device());

    check_cuda(cudaGetLastError());
}

void CudaNonbondedForce::build_lrf_atom_csr(Context& ctx) {
    const int n_group_ranges = data_.group_indices->length;

    n_lrf_source_atom_entries_ = 0;

    lrf_atom_degrees_->zero();
    const int thread_num = 256;
    const int grid_size = (n_group_ranges + thread_num - 1) / thread_num;

    count_lrf_atom_degree_kernel<<<grid_size, thread_num>>>(n_group_ranges, group_pair_modes_->gpu_data_p, data_.group_sizes->gpu_data_p, lrf_atom_degrees_->gpu_data_p);

    check_cuda(cudaGetLastError());

    /*
     * offsets[0] = 0
     * offsets[t+1] = offsets[t] + degrees[t]
     * offsets[n] = total CSR entries
     */

    check_cuda(cub::DeviceScan::ExclusiveSum(
        lrf_scan_temp_->gpu_data_p,
        lrf_scan_temp_bytes_,
        lrf_atom_degrees_->gpu_data_p,
        lrf_atom_offsets_->gpu_data_p,
        n_group_ranges + 1));

    lrf_atom_offsets_->download();
    n_lrf_source_atom_entries_ = lrf_atom_offsets_->cpu_data_p[n_group_ranges];

    if (n_lrf_source_atom_entries_ < 0) {
        throw std::runtime_error("Negative CUDA LRF CSR entry count");
    }

    if (n_lrf_source_atom_entries_ == 0) {
        return;
    }

    const size_t required_capacity = static_cast<size_t>(n_lrf_source_atom_entries_);

    if (!lrf_source_atom_slots_ || required_capacity > lrf_source_atom_capacity_) {
        size_t new_capacity = required_capacity;

        if (lrf_source_atom_capacity_ != 0) {
            const size_t grown_capacity = lrf_source_atom_capacity_ + lrf_source_atom_capacity_ / 2;

            if (grown_capacity > new_capacity) {
                new_capacity = grown_capacity;
            }
        }

        lrf_source_atom_slots_ = std::make_unique<HostDeviceBuffer<int>>(new_capacity, false, true);

        lrf_dx_ = std::make_unique<HostDeviceBuffer<double>>(new_capacity, false, true);

        lrf_dy_ = std::make_unique<HostDeviceBuffer<double>>(new_capacity, false, true);

        lrf_dz_ = std::make_unique<HostDeviceBuffer<double>>(new_capacity, false, true);

        lrf_q_r1_ = std::make_unique<HostDeviceBuffer<double>>(new_capacity, false, true);

        lrf_q_r3_ = std::make_unique<HostDeviceBuffer<double>>(new_capacity, false, true);

        lrf_q_r5_ = std::make_unique<HostDeviceBuffer<double>>(new_capacity, false, true);

        lrf_q_r7_ = std::make_unique<HostDeviceBuffer<double>>(new_capacity, false, true);

        lrf_source_atom_capacity_ = new_capacity;
    }
    fill_lrf_atom_csr_kernel<<<grid_size, thread_num>>>(
        n_group_ranges,
        group_pair_modes_->gpu_data_p,
        data_.group_start_idx->gpu_data_p,
        data_.group_sizes->gpu_data_p,
        lrf_atom_offsets_->gpu_data_p,
        lrf_source_atom_slots_->gpu_data_p);

    check_cuda(cudaGetLastError());
}

void CudaNonbondedForce::build_exact_atom_tiles(Context& ctx) {
    const int n_group_ranges = data_.group_indices->length;

    n_exact_tiles_ = 0;
    n_exact_source_atoms_ = 0;
    if (n_group_ranges == 0) return;

    exact_atom_degrees_->zero();
    exact_atom_offsets_->zero();

    exact_entry_degrees_->zero();
    exact_entry_offsets_->zero();

    const int threads = 256;

    const int blocks = (n_group_ranges + threads - 1) / threads;

    count_exact_atom_tiles_kernel<<<blocks, threads>>>(n_group_ranges,
                                                       group_pair_modes_->gpu_data_p,
                                                       data_.group_start_idx->gpu_data_p,
                                                       data_.group_sizes->gpu_data_p,
                                                       data_.category->gpu_data_p,
                                                       data_.q_state->gpu_data_p,
                                                       exact_atom_degrees_->gpu_data_p,
                                                       exact_entry_degrees_->gpu_data_p);
    check_cuda(cudaGetLastError());

    check_cuda(cub::DeviceScan::ExclusiveSum(
        lrf_scan_temp_->gpu_data_p,
        lrf_scan_temp_bytes_,

        exact_atom_degrees_->gpu_data_p,
        exact_atom_offsets_->gpu_data_p,

        n_group_ranges + 1));

    check_cuda(cub::DeviceScan::ExclusiveSum(
        lrf_scan_temp_->gpu_data_p,
        lrf_scan_temp_bytes_,

        exact_entry_degrees_->gpu_data_p,
        exact_entry_offsets_->gpu_data_p,

        n_group_ranges + 1));

    exact_atom_offsets_->download();
    exact_entry_offsets_->download();

    n_exact_source_atoms_ = exact_atom_offsets_->cpu_data_p[n_group_ranges];

    n_exact_tiles_ = exact_entry_offsets_->cpu_data_p[n_group_ranges];

    if (n_exact_source_atoms_ < 0 || n_exact_tiles_ < 0) {
        throw std::runtime_error("Negative CUDA exact-list size");
    }

    if (static_cast<size_t>(n_exact_tiles_) > exact_tile_capacity_) {
        throw std::runtime_error("CUDA packed exact-tile capacity exceeded");
    }

    const size_t required_atom_capacity = static_cast<size_t>(n_exact_source_atoms_);

    if (required_atom_capacity > exact_source_atom_capacity_) {
        size_t new_capacity = required_atom_capacity;

        if (exact_source_atom_capacity_ != 0) {
            const size_t grown_capacity = exact_source_atom_capacity_ + exact_source_atom_capacity_ / 2;

            new_capacity = std::max(new_capacity, grown_capacity);
        }

        exact_source_atom_slots_ = std::make_unique<HostDeviceBuffer<int>>(new_capacity, false, true);

        exact_source_atom_capacity_ = new_capacity;
    }

    if (n_exact_tiles_ == 0) {
        return;
    }

    fill_exact_atom_tiles_kernel<<<blocks, threads>>>(
        n_group_ranges,

        group_pair_modes_->gpu_data_p,
        data_.group_start_idx->gpu_data_p,
        data_.group_sizes->gpu_data_p,
        data_.category->gpu_data_p,
        data_.q_state->gpu_data_p,

        exact_atom_offsets_->gpu_data_p,
        exact_entry_offsets_->gpu_data_p,

        exact_source_atom_slots_->gpu_data_p,
        exact_tiles_->gpu_data_p);

    check_cuda(cudaGetLastError());
}

void CudaNonbondedForce::calc(Context& ctx) {
    /*
    Sync the coords to CudaNonbondedForce::coords first.
    */
    int sz = data_.n_total;
    int sync_block = 256;
    int sync_grid = (sz + sync_block - 1) / sync_block;
    update_nonbonded_coords_kernel<<<sync_grid, sync_block>>>(ctx.coords->gpu_data_p, data_.atom_idx->gpu_data_p, coord_x_->gpu_data_p, coord_y_->gpu_data_p, coord_z_->gpu_data_p, sz);

    /*
    Do calculation
    */

    if (!ctx.md.lrf || ctx.md.non_bond == 0) {
        calc_all_direct_pairs(ctx);
        return;
    }

    if (ctx.step == ctx.md.steps || ctx.step % ctx.md.non_bond == 0) {
        init_calculation_groups(ctx);
        init_lrf_coefficients(ctx);
    }

    calc_exact_tiles(ctx);
    calc_lrf(ctx);
}