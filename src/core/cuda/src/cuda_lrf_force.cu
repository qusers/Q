#include "cuda_force_accumulation.cuh"
#include "cuda_lrf_force.cuh"
#include "cuda_nonbonded_pair.cuh"
#include "geometry.h"

namespace {

__device__ __forceinline__ real_t combined_slot_distance2(
    int slot1,
    int slot2,

    const real_t* __restrict__ cx,
    const real_t* __restrict__ cy,
    const real_t* __restrict__ cz) {
    const real_t dx = cx[slot1] - cx[slot2];
    const real_t dy = cy[slot1] - cy[slot2];
    const real_t dz = cz[slot1] - cz[slot2];

    return dx * dx + dy * dy + dz * dz;
}

__device__ __forceinline__ real_t combined_slot_center_distance2(
    int slot,
    const coord_t center,

    const real_t* __restrict__ cx,
    const real_t* __restrict__ cy,
    const real_t* __restrict__ cz) {
    const real_t dx = cx[slot] - center.x;
    const real_t dy = cy[slot] - center.y;
    const real_t dz = cz[slot] - center.z;

    return dx * dx + dy * dy + dz * dz;
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

__global__ void build_exact_cluster_tiles_kernel(
    int n_clusters,
    int n_group_ranges,
    int max_exact_tiles,

    const uint8_t* group_pair_modes,
    const int* slot_to_group_range,

    const int* atom_idx,
    const uint8_t* category,
    const int* q_state,
    const uint32_t* pair_masks,

    int* exact_tile_count,
    int* list_overflow,

    ExactEntry* exact_entries,
    uint32_t* exact_pair_masks) {
    constexpr unsigned FULL_MASK = 0xffffffffu;
    constexpr uint8_t Q = static_cast<uint8_t>(AtomCategory::Q);

    const int lane = threadIdx.x & 31;
    const int warp_in_block = threadIdx.x >> 5;
    const int warps_per_block = blockDim.x >> 5;

    const int candidate_index = blockIdx.x * warps_per_block + warp_in_block;

    const int candidate_count = n_clusters * (n_clusters + 1) / 2;

    // Uniform return for the entire warp.
    if (candidate_index >= candidate_count) {
        return;
    }

    const int2 cluster_pair = get_tile_idx(n_clusters, candidate_index);

    const int x_cluster = cluster_pair.x;
    const int y_cluster = cluster_pair.y;

    const int x_start = x_cluster * 32;
    const int y_start = y_cluster * 32;

    const int x_slot = x_start + lane;
    const int x_range = slot_to_group_range[x_slot];
    const int x_atom = atom_idx[x_slot];

    uint32_t pair_mask = 0;

    if (x_range >= 0 && x_atom >= 0) {
        const uint8_t x_category = category[x_slot];
        const int x_state = q_state[x_slot];

#pragma unroll
        for (int y_lane = 0; y_lane < 32; ++y_lane) {
            // Unique ownership within a diagonal tile.
            if (x_cluster == y_cluster && y_lane <= lane) {
                continue;
            }

            const int y_slot = y_start + y_lane;
            const int y_range = slot_to_group_range[y_slot];

            if (y_range < 0) {
                continue;
            }

            const int group_pair_index = get_pair_index(n_group_ranges, x_range, y_range);

            if (group_pair_modes[group_pair_index] != GROUP_PAIR_EXACT) {
                continue;
            }

            const int y_atom = atom_idx[y_slot];

            if (y_atom < 0 || x_atom == y_atom) {
                continue;
            }

            const uint8_t y_category = category[y_slot];

            if (x_category == Q && y_category == Q && x_state != q_state[y_slot]) {
                continue;
            }

            const uint32_t bit = uint32_t{1} << y_lane;
            const size_t row = static_cast<size_t>(candidate_index) * 32 + lane;
            if ((pair_masks[row] & bit) == 0) continue;
            pair_mask |= bit;
        }
    }

    // Test activity after removing exclusions.
    if (!__any_sync(FULL_MASK, pair_mask != 0)) {
        return;
    }

    int output_index = -1;

    if (lane == 0) {
        output_index = atomicAdd(exact_tile_count, 1);

        if (output_index >= max_exact_tiles) {
            atomicExch(list_overflow, 1);
        } else {
            exact_entries[output_index] = {
                x_start, y_start};
        }
    }

    output_index = __shfl_sync(FULL_MASK, output_index, 0);

    if (output_index < max_exact_tiles) {
        const size_t mask_index = static_cast<size_t>(output_index) * 32 + lane;

        exact_pair_masks[mask_index] = pair_mask;
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

__global__ void classify_group_pairs_kernel(
    int n_group_ranges,
    int use_switch_atom,

    real_t solute_solute_cutoff2,
    real_t solute_solvent_cutoff2,
    real_t solvent_solvent_cutoff2,
    real_t rcq2,
    real_t lrf_cutoff2,

    const coord_t solute_center,

    const int* __restrict__ group_start_idx,
    const int* __restrict__ group_sizes,

    const uint8_t* __restrict__ category,
    const int* __restrict__ q_state,

    const real_t* __restrict__ cx,
    const real_t* __restrict__ cy,
    const real_t* __restrict__ cz,

    uint8_t* __restrict__ group_pair_modes) {
    const int pair_index = blockIdx.x * blockDim.x + threadIdx.x;

    const int total_pairs = n_group_ranges * (n_group_ranges + 1) / 2;

    if (pair_index >= total_pairs) {
        return;
    }

    const int2 pair = get_tile_idx(n_group_ranges, pair_index);

    const int range1 = pair.x;
    const int range2 = pair.y;

    const int start1 = group_start_idx[range1];
    const int start2 = group_start_idx[range2];

    const int size1 = group_sizes[range1];
    const int size2 = group_sizes[range2];

    const uint8_t category1 = category[start1];
    const uint8_t category2 = category[start2];

    constexpr uint8_t P = static_cast<uint8_t>(AtomCategory::P);

    constexpr uint8_t Q = static_cast<uint8_t>(AtomCategory::Q);

    constexpr uint8_t W = static_cast<uint8_t>(AtomCategory::W);

    /*
     * Default classification.
     */
    group_pair_modes[pair_index] = GROUP_PAIR_IGNORE;

    /*
     * A group interacts exactly with itself.
     */
    if (range1 == range2) {
        group_pair_modes[pair_index] = GROUP_PAIR_EXACT;
        return;
    }

    const bool range1_is_q = category1 == Q;
    const bool range2_is_q = category2 == Q;

    /*
     * Q-Q is exact only for the same duplicated Q state.
     */
    if (range1_is_q && range2_is_q) {
        if (q_state[start1] == q_state[start2]) {
            group_pair_modes[pair_index] = GROUP_PAIR_EXACT;
        }

        return;
    }

    /*
     * Q-P and Q-W use the rcq rule.
     *
     * With switch atoms:
     *     test the environment switch atom.
     *
     * Without switch atoms:
     *     W still uses its switch atom;
     *     P is exact when any atom is inside rcq.
     */
    if (range1_is_q || range2_is_q) {
        const int environment_range = range1_is_q ? range2 : range1;

        const int environment_start = group_start_idx[environment_range];

        const int environment_size = group_sizes[environment_range];

        const uint8_t environment_category = category[environment_start];

        bool inside_rcq = false;

        if (use_switch_atom != 0 || environment_category == W) {
            inside_rcq = combined_slot_center_distance2(environment_start, solute_center, cx, cy, cz) <= rcq2;
        } else {
            /*
             * The environment is a P group. Check every atom,
             * matching CpuNonbondedForce::inside_rcq().
             */
            for (int local_atom = 0; local_atom < environment_size; ++local_atom) {
                const int slot = environment_start + local_atom;

                if (combined_slot_center_distance2(slot, solute_center, cx, cy, cz) <= rcq2) {
                    inside_rcq = true;
                    break;
                }
            }
        }

        if (inside_rcq) {
            group_pair_modes[pair_index] = GROUP_PAIR_EXACT;
        }

        return;
    }

    /*
     * At this point the pair is P-P, P-W, or W-W.
     */
    real_t normal_cutoff2;

    if (category1 == P && category2 == P) {
        normal_cutoff2 = solute_solute_cutoff2;
    } else if (category1 == W && category2 == W) {
        normal_cutoff2 = solvent_solvent_cutoff2;
    } else {
        normal_cutoff2 = solute_solvent_cutoff2;
    }

    /*
     * Switch-atom mode uses the first slot because the combined
     * list always places the switch atom first in each group.
     */
    if (use_switch_atom != 0) {
        const real_t distance2 = combined_slot_distance2(start1, start2, cx, cy, cz);

        if (distance2 <= normal_cutoff2) {
            group_pair_modes[pair_index] = GROUP_PAIR_EXACT;
        } else if (distance2 <= lrf_cutoff2) {
            group_pair_modes[pair_index] = GROUP_PAIR_LRF;
        }

        return;
    }

    /*
     * Without switch atoms, W-W still uses switch atoms.
     */
    if (category1 == W && category2 == W) {
        const real_t distance2 = combined_slot_distance2(start1, start2, cx, cy, cz);

        if (distance2 <= normal_cutoff2) {
            group_pair_modes[pair_index] = GROUP_PAIR_EXACT;
        } else if (distance2 <= lrf_cutoff2) {
            group_pair_modes[pair_index] = GROUP_PAIR_LRF;
        }

        return;
    }

    /*
     * Without switch atoms, P-W uses the minimum distance from
     * every atom in the P group to the W switch atom.
     */
    if (category1 == W || category2 == W) {
        const int solute_start = category1 == P ? start1 : start2;

        const int solute_size = category1 == P ? size1 : size2;

        const int water_switch = category1 == W ? start1 : start2;

        real_t minimum_distance2 = lrf_cutoff2 + 1.0;

        for (int local_atom = 0; local_atom < solute_size; ++local_atom) {
            const real_t distance2 = combined_slot_distance2(solute_start + local_atom, water_switch, cx, cy, cz);

            /*
             * Exact is the highest-priority classification, so
             * stop as soon as any pair enters the exact cutoff.
             */
            if (distance2 <= normal_cutoff2) {
                group_pair_modes[pair_index] = GROUP_PAIR_EXACT;
                return;
            }

            if (distance2 < minimum_distance2) {
                minimum_distance2 = distance2;
            }
        }

        if (minimum_distance2 <= lrf_cutoff2) {
            group_pair_modes[pair_index] = GROUP_PAIR_LRF;
        }

        return;
    }

    /*
     * P-P without switch atoms:
     * minimum distance over every atom pair.
     */
    real_t minimum_distance2 = lrf_cutoff2 + 1.0;

    for (int local1 = 0; local1 < size1; ++local1) {
        const int slot1 = start1 + local1;

        for (int local2 = 0; local2 < size2; ++local2) {
            const int slot2 = start2 + local2;

            const real_t distance2 = combined_slot_distance2(slot1, slot2, cx, cy, cz);

            if (distance2 <= normal_cutoff2) {
                group_pair_modes[pair_index] = GROUP_PAIR_EXACT;
                return;
            }

            if (distance2 < minimum_distance2) {
                minimum_distance2 = distance2;
            }
        }
    }

    if (minimum_distance2 <= lrf_cutoff2) {
        group_pair_modes[pair_index] = GROUP_PAIR_LRF;
    }
}

__global__ void exact_tiles_nonbonded_force_kernel(
    int n_exact_tiles,
    const ExactEntry* exact_entries,
    const uint32_t* exact_pair_masks,

    int n_states,  // ctx.n_lambdas, used by nb_coul_slot

    // ---- per-atom arrays (length sz, parallel to atom_idx) ----
    const int* atom_idx,               // data_->atom_idx, local i -> global atom index
    const uint8_t* category,           // data_->category, P/Q/W
    const int* q_state,                // data_->q_state, Q state; -1 for P/W
    const real_t* atom_lambdas,        // data_->atom_lambdas
    const real_t* atom_charge,         // data_->atom_charge
    const vdw_atom_param_t* atom_vdw,  // data_->atom_vdw

    // ---- exclusion data ----

    // ---- topology scalars (passed by value) ----
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
    const uint32_t pair_mask = exact_pair_masks[static_cast<size_t>(tile_index) * 32 + lane];

    /*
     * Each lane loads one X atom and one Y atom.
     * Invalid lanes use atom == -1 and still participate in all shuffles.
     */
    const int x_idx = tile.x_start + lane;
    const int y_idx = tile.y_start + lane;

    nonbonded_force_calculation(x_idx, y_idx, pair_mask, tile.x_start, tile.y_start, n_states, atom_idx,
                                category, q_state, atom_lambdas, atom_charge, atom_vdw, coulomb_constant, vdw_rule, cx, cy, cz, dvelocities, e);
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
    real_t sums[LRF_UNIQUE_COMPONENTS];

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

        const real_t charge = atom_charge[source_slot];

        const real_t x = cx[source_slot] - center.x;

        const real_t y = cy[source_slot] - center.y;

        const real_t z = cz[source_slot] - center.z;

        const real_t r2 = x * x + y * y + z * z;

        /*
         * A self relationship should not be LRF, but guard against
         * singular values in case classification changes later.
         */
        if (r2 == 0.0) {
            continue;
        }

        const real_t inv_r = rsqrt(r2);
        const real_t inv_r2 = 1.0 / r2;

        const real_t inv_r3 = inv_r * inv_r2;

        const real_t inv_r5 = inv_r3 * inv_r2;

        const real_t inv_r7 = inv_r5 * inv_r2;

        const real_t q_r1 = charge * inv_r;

        const real_t q_r3 = charge * inv_r3;

        const real_t q_r5 = charge * inv_r5;

        const real_t q_r7 = charge * inv_r7;

        const real_t three_q_r5 = 3.0 * q_r5;

        const real_t nine_q_r5 = 9.0 * q_r5;

        const real_t fifteen_q_r7 = 15.0 * q_r7;

        const real_t xx = x * x;
        const real_t yy = y * y;
        const real_t zz = z * z;

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

void CudaLrfForce::init(Context& ctx, const NonbondedData& data) {
    data_ = &data;
    const int n_group_ranges = data_->group_indices->length;

    const int max_group_pairs = n_group_ranges * (n_group_ranges + 1) / 2;

    if ((data_->n_total & 31) != 0) {
        throw std::runtime_error("CUDA exact cluster list requires n_total to be padded to 32");
    }

    const int n_exact_clusters = data_->n_total / 32;

    const size_t max_exact_tiles = static_cast<size_t>(n_exact_clusters) * static_cast<size_t>(n_exact_clusters + 1) / 2;

    exact_tile_capacity_ = max_exact_tiles;

    group_pair_modes_ = std::make_unique<HostDeviceBuffer<uint8_t>>(max_group_pairs, false, true);

    exact_tiles_ = std::make_unique<HostDeviceBuffer<ExactEntry>>(exact_tile_capacity_, false, true);

    exact_pair_masks_ = std::make_unique<HostDeviceBuffer<uint32_t>>(exact_tile_capacity_ * 32, false, true);

    exact_tile_count_ = std::make_unique<HostDeviceBuffer<int>>(1, true, true);

    list_overflow_ = std::make_unique<HostDeviceBuffer<int>>(1, true, true);

    lrf_coefficients_ = std::make_unique<HostDeviceBuffer<LrfCoefficients>>(ctx.charge_group_config.charge_groups.size(), false, true);

    std::vector<int> slot_to_group_range(data_->n_total, -1);

    for (int range = 0; range < n_group_ranges; ++range) {
        const int begin = data_->group_start_idx->cpu_data_p[range];

        const int size = data_->group_sizes->cpu_data_p[range];

        for (int local_atom = 0; local_atom < size; ++local_atom) {
            slot_to_group_range[begin + local_atom] = range;
        }
    }

    lrf_slot_to_group_range_ = HostDeviceBuffer<int>::from_vector(slot_to_group_range, ctx.command_info.requested_gpu);
}

void CudaLrfForce::init_calculation_groups(Context& ctx, const NonbondedCoordinates& coords) {
    const double solute_solute_cutoff2 = ctx.md.solute_solute * ctx.md.solute_solute;

    const double solute_solvent_cutoff2 = ctx.md.solute_solvent * ctx.md.solute_solvent;

    const double solvent_solvent_cutoff2 = ctx.md.solvent_solvent * ctx.md.solvent_solvent;

    const double rcq2 = ctx.md.q_atom * ctx.md.q_atom;

    const double lrf_cutoff2 = ctx.md.lrf_cutoff * ctx.md.lrf_cutoff;

    const int n_group_ranges = static_cast<int>(data_->group_indices->length);

    if (n_group_ranges == 0) {
        n_exact_tiles_ = 0;
        return;
    }

    constexpr int threads = 256;

    const int total_pairs = n_group_ranges * (n_group_ranges + 1) / 2;

    const int blocks = (total_pairs + threads - 1) / threads;

    int use_switch_atom = ctx.charge_group_config.iuse_switch_atom == 1 ? 1 : 0;

    classify_group_pairs_kernel<<<blocks, threads>>>(
        n_group_ranges,
        use_switch_atom,

        solute_solute_cutoff2,
        solute_solvent_cutoff2,
        solvent_solvent_cutoff2,
        rcq2,
        lrf_cutoff2,

        ctx.topo.solute_center,

        data_->group_start_idx->gpu_data_p,
        data_->group_sizes->gpu_data_p,

        data_->category->gpu_data_p,
        data_->q_state->gpu_data_p,

        coords.x,
        coords.y,
        coords.z,

        group_pair_modes_->gpu_data_p);

    check_cuda(cudaGetLastError());

    build_exact_atom_tiles(ctx);
}

void CudaLrfForce::calc_exact_tiles(Context& ctx, const NonbondedCoordinates& coords) {
    if (n_exact_tiles_ <= 0) {
        return;
    }

    constexpr int thread_num = 256;
    constexpr int warps_per_block = thread_num / 32;

    const int grid_sz = (n_exact_tiles_ + warps_per_block - 1) / warps_per_block;

    exact_tiles_nonbonded_force_kernel<<<grid_sz, thread_num>>>(
        n_exact_tiles_,
        exact_tiles_->gpu_data_p,
        exact_pair_masks_->gpu_data_p,

        ctx.n_lambdas(),

        data_->atom_idx->gpu_data_p,
        data_->category->gpu_data_p,
        data_->q_state->gpu_data_p,
        data_->atom_lambdas->gpu_data_p,
        data_->atom_charge->gpu_data_p,
        data_->atom_vdw->gpu_data_p,

        static_cast<real_t>(ctx.topo.coulomb_constant),

        ctx.topo.vdw_rule,

        coords.x,
        coords.y,
        coords.z,

        ctx.dvelocities->gpu_data_p,
        ctx.energy.device());

    check_cuda(cudaGetLastError());
}

void CudaLrfForce::init_lrf_coefficients(Context& ctx, const NonbondedCoordinates& coords) {
    lrf_coefficients_->zero();

    const int n_group_ranges = static_cast<int>(data_->group_indices->length);

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

        data_->group_indices->gpu_data_p,
        data_->group_start_idx->gpu_data_p,
        data_->group_sizes->gpu_data_p,
        data_->category->gpu_data_p,

        coords.x,
        coords.y,
        coords.z,

        lrf_coefficients_->gpu_data_p);

    check_cuda(cudaGetLastError());
    

    /*
     * Exactly one block per target range.
     */
    build_lrf_coefficients_dense_kernel<<<n_group_ranges, LRF_COEFFICIENT_THREADS>>>(
        n_group_ranges,
        data_->n_total,

        data_->group_indices->gpu_data_p,
        data_->group_start_idx->gpu_data_p,
        data_->category->gpu_data_p,

        group_pair_modes_->gpu_data_p,

        lrf_slot_to_group_range_
            ->gpu_data_p,

        data_->atom_charge->gpu_data_p,

        coords.x,
        coords.y,
        coords.z,

        lrf_coefficients_->gpu_data_p);

    check_cuda(cudaGetLastError());
}

void CudaLrfForce::calc_lrf(Context& ctx, const NonbondedCoordinates& coords) {
    const int n_slots = data_->n_total;

    if (n_slots <= 0) {
        return;
    }

    constexpr int thread_num = 256;

    const int grid_sz = (n_slots + thread_num - 1) / thread_num;

    calc_lrf_kernel<<<grid_sz, thread_num>>>(
        n_slots,

        data_->atom_idx->gpu_data_p,
        data_->atom_to_group->gpu_data_p,
        data_->category->gpu_data_p,
        data_->atom_charge->gpu_data_p,

        coords.x,
        coords.y,
        coords.z,

        lrf_coefficients_->gpu_data_p,

        ctx.topo.coulomb_constant,

        ctx.dvelocities->gpu_data_p,
        ctx.energy.device());

    check_cuda(cudaGetLastError());
}

void CudaLrfForce::build_exact_atom_tiles(Context& ctx) {
    n_exact_tiles_ = 0;

    const int n_group_ranges = static_cast<int>(data_->group_indices->length);

    const int n_clusters = data_->n_total / 32;

    if (n_group_ranges == 0 || n_clusters == 0) {
        return;
    }

    exact_tile_count_->zero();
    list_overflow_->zero();

    constexpr int threads = 256;
    constexpr int warps_per_block = threads / 32;

    const int candidate_count = n_clusters * (n_clusters + 1) / 2;

    const int blocks = (candidate_count + warps_per_block - 1) / warps_per_block;

    build_exact_cluster_tiles_kernel<<<blocks, threads>>>(
        n_clusters,
        n_group_ranges,
        static_cast<int>(exact_tile_capacity_),

        group_pair_modes_->gpu_data_p,
        lrf_slot_to_group_range_->gpu_data_p,

        data_->atom_idx->gpu_data_p,
        data_->category->gpu_data_p,
        data_->q_state->gpu_data_p,
        data_->pair_masks->gpu_data_p,

        exact_tile_count_->gpu_data_p,
        list_overflow_->gpu_data_p,

        exact_tiles_->gpu_data_p,
        exact_pair_masks_->gpu_data_p);

    check_cuda(cudaGetLastError());

    exact_tile_count_->download();
    list_overflow_->download();

    if (list_overflow_->cpu_data_p[0] != 0) {
        throw std::runtime_error("CUDA exact cluster tile capacity exceeded");
    }

    n_exact_tiles_ = exact_tile_count_->cpu_data_p[0];

    if (n_exact_tiles_ < 0 || static_cast<size_t>(n_exact_tiles_) > exact_tile_capacity_) {
        throw std::runtime_error("Invalid CUDA exact cluster tile count");
    }
}

void CudaLrfForce::calc(Context& ctx, const NonbondedCoordinates& coords) {
    if (ctx.step == ctx.md.steps || ctx.step % ctx.md.non_bond == 0) {
        init_calculation_groups(ctx, coords);
        init_lrf_coefficients(ctx, coords);
    }
    calc_exact_tiles(ctx, coords);
    calc_lrf(ctx, coords);
}
