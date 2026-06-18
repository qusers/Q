
#pragma once
#include <cstdint>
#include <memory>

#include "context.h"
#include "host_device_buffer.h"

enum class AtomCategory : uint8_t { P,
                                    Q,
                                    W,
                                    INVALID };
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