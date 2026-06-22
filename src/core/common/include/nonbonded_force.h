
#pragma once
#include <cstdint>
#include <memory>

#include "context.h"
#include "host_device_buffer.h"

// CPU host code and CUDA both include this; keep the slot math in ONE place.
#ifdef __CUDACC__
#define NB_HD __host__ __device__
#else
#define NB_HD
#endif

enum class AtomCategory : uint8_t { P,
                                    Q,
                                    W,
                                    INVALID };

enum NbGroup { NB_PP = 0,
               NB_PW = 1,
               NB_WW = 2 };

NB_HD inline int nb_total_slots(int n_states) { return 3 + 3 * n_states; }
NB_HD inline int nb_qq_slot(int s, int n_states) { return 3 + 0 * n_states + s; }
NB_HD inline int nb_qp_slot(int s, int n_states) { return 3 + 1 * n_states + s; }
NB_HD inline int nb_qw_slot(int s, int n_states) { return 3 + 2 * n_states + s; }

NB_HD inline int nb_energy_slot(uint8_t t1, uint8_t t2,
                                int s1, int s2, int n_states) {
    const uint8_t P = (uint8_t)AtomCategory::P;
    const uint8_t Q = (uint8_t)AtomCategory::Q;
    const uint8_t W = (uint8_t)AtomCategory::W;

    bool q1 = (t1 == Q), q2 = (t2 == Q);
    if (!q1 && !q2) {
        if (t1 == P && t2 == P) return NB_PP;
        if (t1 == W && t2 == W) return NB_WW;
        return NB_PW;  // 其余即 pw
    }
    int state = q1 ? s1 : s2;  // 选 Q 原子的 state
    if (q1 && q2) return nb_qq_slot(state, n_states);
    if (t1 == W || t2 == W) return nb_qw_slot(state, n_states);
    return nb_qp_slot(state, n_states);  // 否则 qp
}

struct NonbondedData {
    int n_total = 0;
    std::unique_ptr<HostDeviceBuffer<int>> atom_idx;      // global atom index
    std::unique_ptr<HostDeviceBuffer<int>> charge_types;  // parallel to atom_idx
    std::unique_ptr<HostDeviceBuffer<int>> catype_types;
    std::unique_ptr<HostDeviceBuffer<uint8_t>> category;     // Atom Category
    std::unique_ptr<HostDeviceBuffer<int>> q_state;          // segment idx; -1 for P/W
    std::unique_ptr<HostDeviceBuffer<real_t>> atom_lambdas;  // lambdas[state]; 1.0 for P/W

    // Precomputed pairwise tables (interned-type indexed) + their descriptors.
    std::unique_ptr<HostDeviceBuffer<real_t>> charge_pair_products;          // [n_ct * n_ct]
    std::unique_ptr<HostDeviceBuffer<vdw_pair_param_t>> catype_pair_params;  // [n_cat * n_cat]

    int n_charge_types = 0;
    int zero_charge_type = -1;
    int n_catype_types = 0;
    int zero_catype_type = -1;

    bool enabled() const { return n_total > 0; }
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
};