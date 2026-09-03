
#pragma once
#include <cstdint>
#include <memory>

#include "constants.h"
#include "context.h"
#include "host_device_buffer.h"

enum class AtomCategory : uint8_t { P,
                                    Q,
                                    W,
                                    INVALID };

// Maps a tile's category/state to the EnergyBuffer *coul* slot.
HD inline int nb_coul_slot(uint8_t t1, uint8_t t2, int s1, int s2, int n_states) {
    const uint8_t P = (uint8_t)AtomCategory::P;
    const uint8_t Q = (uint8_t)AtomCategory::Q;
    const uint8_t W = (uint8_t)AtomCategory::W;
    bool q1 = (t1 == Q), q2 = (t2 == Q);
    if (!q1 && !q2) {
        if (t1 == P && t2 == P) return E_NB_PP_COUL;
        if (t1 == W && t2 == W) return E_NB_WW_COUL;
        return E_NB_PW_COUL;  // remaining non-Q pair is P-W
    }
    int state = q1 ? s1 : s2;  // the Q atom's state
    if (q1 && q2) return EnergyBuffer::eq_index(ENERGY_FIXED_COUNT, state, EQ_NB_QQ_COUL);
    if (t1 == W || t2 == W) return EnergyBuffer::eq_index(ENERGY_FIXED_COUNT, state, EQ_NB_QW_COUL);
    return EnergyBuffer::eq_index(ENERGY_FIXED_COUNT, state, EQ_NB_QP_COUL);  // else Q-P
}

// The vdw slot is always coul+1
HD inline int nb_vdw_slot(uint8_t t1, uint8_t t2, int s1, int s2, int n_states) {
    return nb_coul_slot(t1, t2, s1, s2, n_states) + 1;
}

enum class BondType : uint8_t { Bond23,
                                Bond14,
                                NonBond };

HD inline BondType get_bond_type(int n_atoms_solute, const int* LJ_matrix, int atom1, uint8_t atom1_type, int atom2, uint8_t atom2_type) {
    bool w1 = atom1_type == static_cast<uint8_t>(AtomCategory::W);
    bool w2 = atom2_type == static_cast<uint8_t>(AtomCategory::W);

    if (w1 != w2) {
        return BondType::NonBond;  // solute-water never bonded
    }

    if (w1) {
        // both w
        if ((atom1 - n_atoms_solute) / 3 == (atom2 - n_atoms_solute) / 3) {
            return BondType::Bond23;
        }
        return BondType::NonBond;
    }

    int lj_matrix_value = LJ_matrix[atom1 * n_atoms_solute + atom2];
    if (lj_matrix_value == 3) {
        return BondType::Bond23;
    } else if (lj_matrix_value == 1) {
        return BondType::Bond14;
    }
    return BondType::NonBond;
}

HD inline real_t2 calc_electrostatic(real_t qij, real_t coulomb_constant, real_t inv_dis) {
    if (qij == 0) return {0, 0};
    real_t vel = qij * coulomb_constant * inv_dis;  // k * qi * qj / r
    real_t dvel = -vel * inv_dis;                   /// -k * qi * qj / r2
    return {vel, dvel};
}

HD inline real_t2 calc_vdw(const real_t2& pair, real_t inv_dis) {
    if (pair.x == 0 && pair.y == 0) return {0, 0};
    real_t inv_dis3 = inv_dis * inv_dis * inv_dis;
    real_t inv_dis6 = inv_dis3 * inv_dis3;

    real_t V_a = pair.x * inv_dis6 * inv_dis6;  // precomputed A / r^12
    real_t V_b = pair.y * inv_dis6;             // precomputed B / r^6

    real_t vvdw = V_a - V_b;
    real_t dvvdw = (static_cast<real_t>(-12.0) * V_a + static_cast<real_t>(6.0) * V_b) * inv_dis;

    return {vvdw, dvvdw};
}

HD inline real_t2 combine_vdw(int vdw_rule, real_t aii_i, real_t bii_i, real_t aii_j, real_t bii_j) {
    real_t a, b;
    if (vdw_rule == VDW_GEOMETRIC) {
        calc_vdw_geometric(aii_i, aii_j, bii_i, bii_j, (real_t)1.0, &a, &b);
    } else {
        calc_vdw_arithmetic(aii_i, aii_j, bii_i, bii_j, (real_t)1.0, &a, &b);
    }
    return {a, b};
}

struct NonbondedData {
    int n_total = 0;
    std::unique_ptr<HostDeviceBuffer<int>> atom_idx;       // global atom index
    std::unique_ptr<HostDeviceBuffer<int>> atom_to_group;  // global atom index

    std::unique_ptr<HostDeviceBuffer<uint8_t>> category;     // Atom Category
    std::unique_ptr<HostDeviceBuffer<int>> q_state;          // segment idx; -1 for P/W
    std::unique_ptr<HostDeviceBuffer<real_t>> atom_lambdas;  // lambdas[state]; 1.0 for P/W
    std::unique_ptr<HostDeviceBuffer<real_t>> atom_charge;
    std::unique_ptr<HostDeviceBuffer<vdw_atom_param_t>> atom_vdw;

    bool enabled() const { return n_total > 0; }
};

struct LrfCoefficients {
    coord_t center{};

    double phi0 = 0;
    double phi1[3]{};
    double phi2[9]{};
    double phi3[27]{};
};

class NonbondedForce {
   public:
    virtual ~NonbondedForce() = default;
    // Builds the combined atom list + pairwise tables, then calls init_backend(ctx).
    // Call once before the first MD step
    void init(Context& ctx);

    // Per-step : compute all groups, write energies to ctx.E_nonbond_*, ctx.EQ_nonbond_*
    virtual void calc(Context& ctx) = 0;

    virtual void cleanup() {}

   protected:
    virtual void init_backend(Context& ctx) {}  // backend allocation hook

    NonbondedData data_;

   private:
    void build_combinded_list(Context& ctx);  // atom_idx, category, q_state, atom_lambdas
    void build_charge_table(Context& ctx);    // charge_pair_products + charge_types + counts
    void build_catype_table(Context& ctx);    // catype_pair_params + catype_types + counts
    void build_atom_to_group(Context& ctx);   // atom_to_group
};
