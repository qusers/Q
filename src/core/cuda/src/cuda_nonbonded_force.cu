#include <cub/device/device_scan.cuh>

#include "cuda_force_accumulation.cuh"
#include "cuda_nonbonded_force.cuh"
#include "geometry.h"

namespace {

/*
 * Unique component layout:
 *
 *  0  phi0
 *
 *  1  phi1.x
 *  2  phi1.y
 *  3  phi1.z
 *
 *  4  phi2.xxx -> xx
 *  5  phi2.yy
 *  6  phi2.zz
 *  7  phi2.xy
 *  8  phi2.xz
 *  9  phi2.yz
 *
 * 10  phi3.xxx
 * 11  phi3.yyy
 * 12  phi3.zzz
 * 13  phi3.xxy
 * 14  phi3.xxz
 * 15  phi3.xyy
 * 16  phi3.yyz
 * 17  phi3.xzz
 * 18  phi3.yzz
 * 19  phi3.xyz
 */
__device__ __forceinline__ double
calculate_unique_lrf_component(
    int phi,
    double charge,
    double x,
    double y,
    double z) {
    const double r2 =
        x * x + y * y + z * z;

    if (r2 == 0.0) {
        return 0.0;
    }

    const double inv_r = rsqrt(r2);
    const double inv_r2 = 1.0 / r2;
    const double inv_r3 = inv_r * inv_r2;
    const double inv_r5 = inv_r3 * inv_r2;
    const double inv_r7 = inv_r5 * inv_r2;

    switch (phi) {
        case 0:
            return charge * inv_r;

        case 1:
            return -charge * x * inv_r3;

        case 2:
            return -charge * y * inv_r3;

        case 3:
            return -charge * z * inv_r3;

        case 4:
            return charge *
                   (3.0 * x * x * inv_r5 -
                    inv_r3);

        case 5:
            return charge *
                   (3.0 * y * y * inv_r5 -
                    inv_r3);

        case 6:
            return charge *
                   (3.0 * z * z * inv_r5 -
                    inv_r3);

        case 7:
            return charge *
                   (3.0 * x * y * inv_r5);

        case 8:
            return charge *
                   (3.0 * x * z * inv_r5);

        case 9:
            return charge *
                   (3.0 * y * z * inv_r5);

        case 10:
            return charge *
                   (9.0 * x * inv_r5 -
                    15.0 * x * x * x * inv_r7);

        case 11:
            return charge *
                   (9.0 * y * inv_r5 -
                    15.0 * y * y * y * inv_r7);

        case 12:
            return charge *
                   (9.0 * z * inv_r5 -
                    15.0 * z * z * z * inv_r7);

        case 13:
            return charge *
                   (3.0 * y * inv_r5 -
                    15.0 * x * x * y * inv_r7);

        case 14:
            return charge *
                   (3.0 * z * inv_r5 -
                    15.0 * x * x * z * inv_r7);

        case 15:
            return charge *
                   (3.0 * x * inv_r5 -
                    15.0 * x * y * y * inv_r7);

        case 16:
            return charge *
                   (3.0 * z * inv_r5 -
                    15.0 * y * y * z * inv_r7);

        case 17:
            return charge *
                   (3.0 * x * inv_r5 -
                    15.0 * x * z * z * inv_r7);

        case 18:
            return charge *
                   (3.0 * y * inv_r5 -
                    15.0 * y * z * z * inv_r7);

        case 19:
            return charge *
                   (-15.0 * x * y * z * inv_r7);

        default:
            return 0.0;
    }
}

__device__ __forceinline__ void
write_unique_lrf_component(
    LrfCoefficients& output,
    int phi,
    double value) {
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

__device__ void accumulate_lrf_direction(
    int source_range,
    int target_group,

    const int* group_start_idx,
    const int* group_sizes,

    const real_t* atom_charge,
    const real_t* cx,
    const real_t* cy,
    const real_t* cz,

    LrfCoefficients* coefficients) {
    constexpr unsigned FULL_MASK = 0xffffffffu;

    const int lane = threadIdx.x & 31;

    const int source_start = group_start_idx[source_range];

    const int source_size = group_sizes[source_range];

    const coord_t target_center = coefficients[target_group].center;

    double local_phi0 = 0.0;
    double local_phi1[3] = {};
    double local_phi2[9] = {};
    double local_phi3[27] = {};

    for (int local_atom = lane; local_atom < source_size; local_atom += 32) {
        const int slot = source_start + local_atom;

        const double charge = static_cast<double>(atom_charge[slot]);

        const double rx = static_cast<double>(cx[slot]) - target_center.x;

        const double ry = static_cast<double>(cy[slot]) - target_center.y;

        const double rz = static_cast<double>(cz[slot]) - target_center.z;

        const double r[3] = {rx, ry, rz};

        const double r2 = rx * rx + ry * ry + rz * rz;

        const double r_length = sqrt(r2);

        const double inv_r = 1.0 / r_length;

        const double inv_r2 = 1.0 / r2;

        const double inv_r3 = inv_r * inv_r2;

        const double inv_r5 = inv_r3 * inv_r2;

        const double inv_r7 = inv_r5 * inv_r2;

        /*
         * phi0 += q/r
         */
        local_phi0 += charge * inv_r;

        /*
         * phi1[a] -= q*r[a]/r^3
         */
        for (int a = 0; a < 3; ++a) {
            local_phi1[a] -= charge * r[a] * inv_r3;
        }
        /*
         * phi2[a,b] += q *
         *     (3*r[a]*r[b]/r^5 - delta[a,b]/r^3)
         */
        for (int a = 0; a < 3; ++a) {
            for (int b = 0; b < 3; ++b) {
                const int index = a * 3 + b;

                const double delta_ab = a == b ? 1.0 : 0.0;

                local_phi2[index] += charge * (3.0 * r[a] * r[b] * inv_r5 - delta_ab * inv_r3);
            }
        }

        /*
         * phi3[a,b,c] += q * (
         *     3*(delta_ab*r[c] +
         *        delta_ac*r[b] +
         *        delta_bc*r[a])/r^5
         *     - 15*r[a]*r[b]*r[c]/r^7
         * )
         */
        for (int a = 0; a < 3; ++a) {
            for (int b = 0; b < 3; ++b) {
                for (int c = 0; c < 3; ++c) {
                    const int index = (a * 3 + b) * 3 + c;

                    const double delta_ab = a == b ? 1.0 : 0.0;
                    const double delta_ac = a == c ? 1.0 : 0.0;
                    const double delta_bc = b == c ? 1.0 : 0.0;

                    const double v1 = 3.0 * (delta_ab * r[c] + delta_ac * r[b] + delta_bc * r[a]) * inv_r5;

                    const double v2 = -15.0 * r[a] * r[b] * r[c] * inv_r7;

                    local_phi3[index] += charge * (v1 + v2);
                }
            }
        }
    }

    /*
     * Reduce all lane-local coefficients to lane 0.
     */
    for (int offset = 16; offset > 0; offset >>= 1) {
        local_phi0 += __shfl_down_sync(FULL_MASK, local_phi0, offset);

        for (int i = 0; i < 3; ++i) {
            local_phi1[i] += __shfl_down_sync(FULL_MASK, local_phi1[i], offset);
        }

        for (int i = 0; i < 9; ++i) {
            local_phi2[i] += __shfl_down_sync(FULL_MASK, local_phi2[i], offset);
        }

        for (int i = 0; i < 27; ++i) {
            local_phi3[i] += __shfl_down_sync(FULL_MASK, local_phi3[i], offset);
        }
    }

    if (lane == 0) {
        atomicAdd(&coefficients[target_group].phi0, local_phi0);

        for (int i = 0; i < 3; ++i) {
            atomicAdd(&coefficients[target_group].phi1[i], local_phi1[i]);
        }

        for (int i = 0; i < 9; ++i) {
            atomicAdd(&coefficients[target_group].phi2[i], local_phi2[i]);
        }

        for (int i = 0; i < 27; ++i) {
            atomicAdd(&coefficients[target_group].phi3[i], local_phi3[i]);
        }
    }
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

__global__ void build_lrf_coefficients_csr_kernel(
    int n_group_ranges,
    const int* group_indices,
    const int* group_start_idx,
    const uint8_t* category,

    const int* atom_offsets,
    const int* source_atom_slots,

    const real_t* atom_charge,

    const real_t* cx,
    const real_t* cy,
    const real_t* cz,
    LrfCoefficients* coefficients) {
    const int target_range = blockIdx.x;
    const int phi = blockIdx.y;
    const int thread = threadIdx.x;

    if (target_range >= n_group_ranges || phi >= 20) {
        return;
    }

    const int target_start = group_start_idx[target_range];

    const uint8_t target_category = category[target_start];

    constexpr uint8_t P = static_cast<uint8_t>(AtomCategory::P);

    constexpr uint8_t W = static_cast<uint8_t>(AtomCategory::W);

    /*
     * LRF coefficients are only required for non-Q groups.
     * This condition is block-uniform.
     */
    if (target_category != P && target_category != W) {
        return;
    }

    const int target_group = group_indices[target_range];

    const coord_t target_center = coefficients[target_group].center;

    const int begin = atom_offsets[target_range];

    const int end = atom_offsets[target_range + 1];

    double local_sum = 0.0;

    for (int entry = begin + thread; entry < end; entry += blockDim.x) {
        const int slot = source_atom_slots[entry];

        const double charge = static_cast<double>(atom_charge[slot]);

        const double x = static_cast<double>(cx[slot]) - target_center.x;

        const double y = static_cast<double>(cy[slot]) - target_center.y;

        const double z = static_cast<double>(cz[slot]) - target_center.z;

        local_sum += calculate_unique_lrf_component(phi, charge, x, y, z);
    }

    constexpr unsigned FULL_MASK = 0xffffffffu;

    for (int offset = 16; offset > 0; offset >>= 1) {
        local_sum += __shfl_down_sync(FULL_MASK, local_sum, offset);
    }

    __shared__ double warp_sums[4];

    const int lane = thread & 31;
    const int warp = thread >> 5;

    if (lane == 0) {
        warp_sums[warp] = local_sum;
    }

    __syncthreads();

    if (warp == 0) {
        double block_sum = lane < 4 ? warp_sums[lane] : 0.0;

        for (int offset = 16; offset > 0; offset >>= 1) {
            block_sum += __shfl_down_sync(FULL_MASK, block_sum, offset);
        }

        if (lane == 0) {
            write_unique_lrf_component(coefficients[target_group], phi, block_sum);
        }
    }
}

__global__ void build_lrf_coefficients_kernel(
    int n_group_ranges,

    const uint8_t* group_pair_modes,
    const int* group_indices,
    const int* group_start_idx,
    const int* group_sizes,

    const uint8_t* category,
    const real_t* atom_charge,

    const real_t* cx,
    const real_t* cy,
    const real_t* cz,

    LrfCoefficients* coefficients) {
    constexpr unsigned FULL_MASK = 0xffffffffu;

    const int lane = threadIdx.x & 31;
    const int warp_in_block = threadIdx.x >> 5;
    const int warps_per_block = blockDim.x >> 5;
    const int target_range = blockIdx.x * warps_per_block + warp_in_block;

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
    const int target_group = group_indices[target_range];
    const coord_t target_center = coefficients[target_group].center;

    double local_phi0 = 0.0;
    double local_phi1[3] = {};
    double local_phi2[9] = {};
    double local_phi3[27] = {};

    for (int source_range = 0; source_range < n_group_ranges; ++source_range) {
        const int pair_index = get_pair_index(n_group_ranges, target_range, source_range);
        if (group_pair_modes[pair_index] != GROUP_PAIR_LRF) {
            continue;
        }

        const int source_start = group_start_idx[source_range];
        const int source_size = group_sizes[source_range];

        for (int local_atom = lane; local_atom < source_size; local_atom += 32) {
            const int slot = source_start + local_atom;

            const double charge = static_cast<double>(atom_charge[slot]);

            const double rx = static_cast<double>(cx[slot]) - target_center.x;
            const double ry = static_cast<double>(cy[slot]) - target_center.y;
            const double rz = static_cast<double>(cz[slot]) - target_center.z;

            const double r[3] = {rx, ry, rz};

            const double r2 = rx * rx + ry * ry + rz * rz;

            /*
             * A valid LRF pair should never contain the target group
             * itself, but guard against singular input anyway.
             */
            if (r2 == 0.0) {
                continue;
            }

            const double inv_r = rsqrt(r2);
            const double inv_r2 = 1.0 / r2;
            const double inv_r3 = inv_r * inv_r2;
            const double inv_r5 = inv_r3 * inv_r2;
            const double inv_r7 = inv_r5 * inv_r2;

            /*
             * phi0 = sum(q/r)
             */
            local_phi0 += charge * inv_r;

            /*
             * phi1[a] = sum(-q*r[a]/r^3)
             */
#pragma unroll
            for (int a = 0; a < 3; ++a) {
                local_phi1[a] -= charge * r[a] * inv_r3;
            }

            /*
             * phi2[a,b] =
             * q * (3*r[a]*r[b]/r^5 - delta[a,b]/r^3)
             */
#pragma unroll
            for (int a = 0; a < 3; ++a) {
#pragma unroll
                for (int b = 0; b < 3; ++b) {
                    const int index = a * 3 + b;
                    const double delta_ab = a == b ? 1.0 : 0.0;

                    local_phi2[index] += charge * (3.0 * r[a] * r[b] * inv_r5 - delta_ab * inv_r3);
                }
            }

            /*
             * phi3[a,b,c] =
             * q * (
             *   3*(delta_ab*r[c] +
             *      delta_ac*r[b] +
             *      delta_bc*r[a])/r^5
             *   - 15*r[a]*r[b]*r[c]/r^7
             * )
             */
#pragma unroll
            for (int a = 0; a < 3; ++a) {
#pragma unroll
                for (int b = 0; b < 3; ++b) {
#pragma unroll
                    for (int c = 0; c < 3; ++c) {
                        const int index = (a * 3 + b) * 3 + c;

                        const double delta_ab = a == b ? 1.0 : 0.0;
                        const double delta_ac = a == c ? 1.0 : 0.0;
                        const double delta_bc = b == c ? 1.0 : 0.0;

                        const double v1 = 3.0 * (delta_ab * r[c] + delta_ac * r[b] + delta_bc * r[a]) * inv_r5;

                        const double v2 = -15.0 * r[a] * r[b] * r[c] * inv_r7;

                        local_phi3[index] += charge * (v1 + v2);
                    }
                }
            }
        }
    }

    for (int offset = 16; offset > 0; offset >>= 1) {
        local_phi0 += __shfl_down_sync(FULL_MASK, local_phi0, offset);

#pragma unroll
        for (int i = 0; i < 3; ++i) {
            local_phi1[i] += __shfl_down_sync(FULL_MASK, local_phi1[i], offset);
        }

#pragma unroll
        for (int i = 0; i < 9; ++i) {
            local_phi2[i] += __shfl_down_sync(FULL_MASK, local_phi2[i], offset);
        }

#pragma unroll
        for (int i = 0; i < 27; ++i) {
            local_phi3[i] += __shfl_down_sync(FULL_MASK, local_phi3[i], offset);
        }
    }

    if (lane == 0) {
        LrfCoefficients& output = coefficients[target_group];

        /*
         * Do not overwrite output.center, which was initialized by
         * compute_lrf_centers_kernel().
         */
        output.phi0 = local_phi0;

#pragma unroll
        for (int i = 0; i < 3; ++i) {
            output.phi1[i] = local_phi1[i];
        }

#pragma unroll
        for (int i = 0; i < 9; ++i) {
            output.phi2[i] = local_phi2[i];
        }

#pragma unroll
        for (int i = 0; i < 27; ++i) {
            output.phi3[i] = local_phi3[i];
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

__global__ void build_pair_lists_kernel(
    int n_groups_ranges,

    const uint8_t* group_pair_modes,
    const int* group_indices,
    const int* group_start_idx,
    const int* group_sizes,

    const int exact_tile_capacity,
    int* exact_tile_count,
    ExactEntry* exact_tiles,

    const int lrf_pair_capacity,
    int* lrf_pair_count,
    LrfPairEntry* lrf_pairs,

    int* overflow

) {
    const int pair_index = blockIdx.x * blockDim.x + threadIdx.x;
    const int total_pairs = n_groups_ranges * (n_groups_ranges + 1) / 2;
    if (pair_index >= total_pairs) {
        return;
    }
    const int2 pair = get_tile_idx(n_groups_ranges, pair_index);
    const int group1 = pair.x;
    const int group2 = pair.y;

    if (group2 < group1) {
        return;
    }

    const uint8_t mode = group_pair_modes[pair_index];
    if (mode == GROUP_PAIR_IGNORE) {
        return;
    }

    if (mode == GROUP_PAIR_LRF) {
        const int dst = atomicAdd(lrf_pair_count, 1);

        if (dst >= lrf_pair_capacity) {
            atomicExch(overflow, 1);
            return;
        }

        lrf_pairs[dst] = {
            group1,
            group2,
            group_indices[group1],
            group_indices[group2],
        };
        return;
    }

    const int start1 = group_start_idx[group1];
    const int start2 = group_start_idx[group2];
    const int size1 = group_sizes[group1];
    const int size2 = group_sizes[group2];

    const int nx = (size1 + 31) / 32;
    const int ny = (size2 + 31) / 32;

    int tile_count = 0;
    if (group1 == group2) {
        tile_count = nx * (nx + 1) / 2;
    } else {
        tile_count = nx * ny;
    }

    const int base = atomicAdd(exact_tile_count, tile_count);
    if (base + tile_count > exact_tile_capacity) {
        atomicExch(overflow, 1);
        return;
    }

    int dst = base;
    for (int ix = 0; ix < nx; ix++) {
        const int x_offset = ix * 32;
        const int x_len = min(32, size1 - x_offset);

        for (int iy = 0; iy < ny; iy++) {
            if (group1 == group2 && iy < ix) {
                continue;
            }

            const int y_offset = iy * 32;
            const int y_len = min(32, size2 - y_offset);

            ExactEntry tile;
            tile.x_start = start1 + x_offset;
            tile.y_start = start2 + y_offset;
            tile.x_len = x_len;
            tile.y_len = y_len;
            tile.diagonal = (group1 == group2 && ix == iy);
            exact_tiles[dst++] = tile;
        }
    }
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
    const int y_idx = lane < tile.y_len ? tile.y_start + lane : -1;

    nonbonded_force_calculation(x_idx, y_idx, tile.diagonal, tile.x_start, tile.y_start, n_states, n_atoms_solute, atom_idx,
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
    build_pair_lists_kernel<<<grid, thread_num>>>(n_group_ranges,
                                                  group_pair_modes_->gpu_data_p,
                                                  data_.group_indices->gpu_data_p,
                                                  data_.group_start_idx->gpu_data_p,
                                                  data_.group_sizes->gpu_data_p,
                                                  exact_tile_capacity_,
                                                  exact_tile_count_->gpu_data_p,
                                                  exact_tiles_->gpu_data_p,
                                                  lrf_pair_capacity_,
                                                  lrf_pair_count_->gpu_data_p,
                                                  lrf_group_pairs_->gpu_data_p,
                                                  list_overflow_->gpu_data_p);

    check_cuda(cudaGetLastError());
    exact_tile_count_->download();
    lrf_pair_count_->download();
    list_overflow_->download();

    if (list_overflow_->cpu_data_p[0] != 0) {
        throw std::runtime_error("CUDA nonbonded pair-list capacity exceeded");
    }

    n_exact_tiles_ = exact_tile_count_->cpu_data_p[0];
    n_lrf_pairs_ = lrf_pair_count_->cpu_data_p[0];

    if (n_exact_tiles_ < 0 || static_cast<size_t>(n_exact_tiles_) > exact_tile_capacity_) {
        throw std::runtime_error("Invalid CUDA exact tile count");
    }

    if (n_lrf_pairs_ < 0 || static_cast<size_t>(n_lrf_pairs_) > lrf_pair_capacity_) {
        throw std::runtime_error("Invalid CUDA LRF pair count");
    }
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
    build_lrf_atom_csr(ctx);
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
    constexpr int thread_num = 256;
    constexpr int warps_per_block = thread_num / 32;
    const int n_group_ranges = static_cast<int>(data_.group_indices->length);
    const int grid = (n_group_ranges + warps_per_block - 1) / warps_per_block;

    if (n_group_ranges > 0) {
        compute_lrf_centers_kernel<<<grid, thread_num>>>(
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
    }
    if (n_lrf_source_atom_entries_ > 0) {
        constexpr int coefficient_threads = 128;
        constexpr int unique_phi_count = 20;

        dim3 coefficient_grid = (n_group_ranges, unique_phi_count, 1);

        build_lrf_coefficients_csr_kernel<<<coefficient_grid, coefficient_threads>>>(
            n_group_ranges,

            data_.group_indices->gpu_data_p,
            data_.group_start_idx->gpu_data_p,
            data_.category->gpu_data_p,

            lrf_atom_offsets_->gpu_data_p,
            lrf_source_atom_slots_->gpu_data_p,

            data_.atom_charge->gpu_data_p,

            coord_x_->gpu_data_p,
            coord_y_->gpu_data_p,
            coord_z_->gpu_data_p,

            lrf_coefficients_->gpu_data_p);

        check_cuda(cudaGetLastError());
    }
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