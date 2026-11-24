#pragma once
#include "cuda_force_accumulation.cuh"
#include "nonbonded_force.h"

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

__device__ void compute_pair(
    // atom1
    real_t atom1_charge, const vdw_atom_param_t& atom1_vdw, real_t atom1_lambda, const real_t3& atom1_coord,
    /// atom2
    real_t atom2_charge, const vdw_atom_param_t& atom2_vdw, real_t atom2_lambda, const real_t3& atom2_coord,
    // scalars
    real_t coulomb_constant, int vdw_rule,
    // soft-core: soft_tile is false when neither tile holds Q atoms
    bool soft_tile, const softcore_pair_params_t& sc_params, real_t atom1_sc, real_t atom2_sc,
    // foreign work: tile_state is the FEP state of the tile's Q atoms
    const foreign_targets_t& foreign, int tile_state, int n_states, energy_accum_t* e,
    // REST2: pair_rest scales the force of this tile's pairs; q_pair = both atoms are Q atoms
    real_t pair_rest, bool q_pair,
    // output
    real_t3& atom1_force, real_t3& atom2_force,
    real_t& e_coul, real_t& e_vdw) {
    real_t dx = atom2_coord.x - atom1_coord.x;
    real_t dy = atom2_coord.y - atom1_coord.y;
    real_t dz = atom2_coord.z - atom1_coord.z;
    real_t dis2 = dx * dx + dy * dy + dz * dz;
    real_t inv_dis2 = static_cast<real_t>(1.0) / dis2;
    real_t inv_dis = sqrt(inv_dis2);

    real_t qij = atom1_charge * atom2_charge;
    real_t2 pair = combine_vdw(vdw_rule, atom1_vdw.aii_normal, atom1_vdw.bii_normal, atom2_vdw.aii_normal, atom2_vdw.bii_normal);

    real_t2 el = calc_electrostatic(qij, coulomb_constant, inv_dis);
    real_t2 vdw = calc_vdw(pair, inv_dis);

    if (soft_tile) {
        const real_t sc = max(atom1_sc, atom2_sc);  // < 0: neither atom has soft-core
        if (sc >= 0) {
            const real_t sigma = softcore_pair_sigma(vdw_rule, atom1_vdw.aii_normal, atom1_vdw.bii_normal, atom2_vdw.aii_normal, atom2_vdw.bii_normal);
            const real_t dis = sqrt(dis2);
            if (sc > 0) apply_gapsys(sc_params, sc, sigma, qij, coulomb_constant, pair, dis, el, vdw);
            // At r = 0 a soft pair has a finite energy and no force direction.
            if (!(dis2 > 0)) inv_dis = 0;

            for (int k = 0; k < foreign.n; k++) {
                const real_t target_sc = static_cast<real_t>(foreign.sc[k][tile_state]);
                const real_t diff = gapsys_shape_difference(sc_params, sc, target_sc, sigma, qij, coulomb_constant, pair, dis);
                if (diff != 0) {
                    const double target_rest = q_pair ? foreign.rest_target[k] : foreign.rest_target_sqrt[k];
                    atomic_add_energy(&e[EnergyBuffer::foreign_index(n_states, k)], static_cast<real_t>(foreign.weight[k][tile_state] * target_rest) * diff);
                }
            }
        }
    }

    real_t lambda = min(atom1_lambda, atom2_lambda) * pair_rest;

    real_t dva = (el.y + vdw.y) * inv_dis * lambda;

    atom1_force.x -= dva * dx;
    atom1_force.y -= dva * dy;
    atom1_force.z -= dva * dz;

    atom2_force.x += dva * dx;
    atom2_force.y += dva * dy;
    atom2_force.z += dva * dz;

    e_coul += el.x;
    e_vdw += vdw.x;
}

__device__ void shuffle(int& atom, real_t& atom_charge, vdw_atom_param_t& atom_vdw, real_t& atom_lambda, real_t3& atom_force, real_t3& atom_coord, real_t& atom_sc, bool soft_tile) {
    constexpr unsigned FULL_MASK = 0xFFFFFFFF;
    int src = ((threadIdx.x & 31) + 1) & 31;
    atom = __shfl_sync(FULL_MASK, atom, src);
    atom_charge = __shfl_sync(FULL_MASK, atom_charge, src);

    atom_vdw.aii_normal = __shfl_sync(FULL_MASK, atom_vdw.aii_normal, src);
    atom_vdw.bii_normal = __shfl_sync(FULL_MASK, atom_vdw.bii_normal, src);

    atom_lambda = __shfl_sync(FULL_MASK, atom_lambda, src);
    atom_force.x = __shfl_sync(FULL_MASK, atom_force.x, src);
    atom_force.y = __shfl_sync(FULL_MASK, atom_force.y, src);
    atom_force.z = __shfl_sync(FULL_MASK, atom_force.z, src);
    atom_coord.x = __shfl_sync(FULL_MASK, atom_coord.x, src);
    atom_coord.y = __shfl_sync(FULL_MASK, atom_coord.y, src);
    atom_coord.z = __shfl_sync(FULL_MASK, atom_coord.z, src);
    if (soft_tile) atom_sc = __shfl_sync(FULL_MASK, atom_sc, src);
}

// kSoft = false compiles the soft-core and foreign-work code out of the kernel, so a run
// without soft-core pays nothing for it.
template <bool kSoft>
__device__ void nonbonded_force_calculation(
    int x_idx,
    int y_idx,
    uint32_t pair_mask,
    int base_x,
    int base_y,

    int n_states,  // ctx.n_lambdas, used by nb_coul_slot

    // ---- per-atom arrays (length sz, parallel to atom_idx) ----
    const int* atom_idx,               // data_.atom_idx, local i -> global atom index
    const uint8_t* category,           // data_.category, P/Q/W
    const int* q_state,                // data_.q_state, Q state; -1 for P/W
    const real_t* atom_lambdas,        // data_.atom_lambdas
    const real_t* atom_charge,         // data_.atom_charge
    const vdw_atom_param_t* atom_vdw,  // data_.atom_vdw

    // ---- topology scalars (passed by value) ----
    real_t coulomb_constant,  // ctx.topo.coulomb_constant
    int vdw_rule,
    const real_t* atom_softcore,
    softcore_pair_params_t sc_params,
    foreign_targets_t foreign,

    // ---- coordinates / outputs ----
    const real_t* cx, const real_t* cy, const real_t* cz,
    dvel_t* dvelocities,  // ctx.dvelocities->gpu_data_p (fixed-point, atomic_add_force)

    // ---- energy accumulators  ----
    energy_accum_t* e

) {
    int lane = threadIdx.x & 31;
    const int atom1 = x_idx == -1 ? -1 : atom_idx[x_idx];

    const real_t atom1_charge = atom1 == -1 ? 0 : atom_charge[x_idx];
    const vdw_atom_param_t atom1_vdw = atom1 == -1 ? vdw_atom_param_t{} : atom_vdw[x_idx];
    const real_t atom1_lambda = atom1 == -1 ? 0 : atom_lambdas[x_idx];
    real_t3 atom1_coord = atom1 == -1 ? real_t3{0, 0, 0} : real_t3{cx[x_idx], cy[x_idx], cz[x_idx]};
    real_t3 atom1_force = {0, 0, 0};

    int atom2 = y_idx == -1 ? -1 : atom_idx[y_idx];
    real_t atom2_charge = atom2 == -1 ? 0 : atom_charge[y_idx];
    vdw_atom_param_t atom2_vdw = atom2 == -1 ? vdw_atom_param_t{} : atom_vdw[y_idx];
    real_t atom2_lambda = atom2 == -1 ? 0 : atom_lambdas[y_idx];
    real_t3 atom2_coord = atom2 == -1 ? real_t3{0, 0, 0} : real_t3{cx[y_idx], cy[y_idx], cz[y_idx]};
    real_t3 atom2_force = {0, 0, 0};

    // A cluster holds atoms of one category, so the first slot tells whether a tile has Q atoms.
    const uint8_t q_category = static_cast<uint8_t>(AtomCategory::Q);
    const bool q_x = category[base_x] == q_category;
    const bool q_y = category[base_y] == q_category;
    const bool soft_tile = kSoft && (q_x || q_y);
    const bool q_pair = q_x && q_y;
    const real_t pair_rest = static_cast<real_t>(q_pair ? foreign.rest : (q_x || q_y) ? foreign.rest_sqrt
                                                                                      : 1.0);
    const real_t atom1_sc = (soft_tile && atom1 != -1) ? atom_softcore[x_idx] : -1;
    real_t atom2_sc = (soft_tile && atom2 != -1) ? atom_softcore[y_idx] : -1;
    const int tile_state = category[base_x] == q_category ? q_state[base_x] : q_state[base_y];

    real_t local_e_coul = 0, local_e_vdw = 0;
    int y_lane = lane;
    for (int i = 0; i < 32; i++) {
        const uint32_t bit = uint32_t{1} << y_lane;
        const bool mask_enabled = (pair_mask & bit) != 0;
        if (mask_enabled) {
            compute_pair(atom1_charge, atom1_vdw, atom1_lambda, atom1_coord,
                         atom2_charge, atom2_vdw, atom2_lambda, atom2_coord,
                         coulomb_constant, vdw_rule,
                         soft_tile, sc_params, atom1_sc, atom2_sc,
                         foreign, tile_state, n_states, e,
                         pair_rest, q_pair,
                         atom1_force, atom2_force,
                         local_e_coul, local_e_vdw);
        }
        shuffle(atom2, atom2_charge, atom2_vdw, atom2_lambda, atom2_force, atom2_coord, atom2_sc, soft_tile);
        y_lane = (y_lane + 1) & 31;
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

}  // namespace
