#pragma once
#include <memory>
#include <set>

#include "context.h"
#include "host_device_buffer.h"
#include "md_types.h"

struct ConstraintBond {
    int ai;
    int aj;
    double dist2;
};

struct ConstraintData {
    int n_constraints = 0;
    std::unique_ptr<HostDeviceBuffer<int>> mol_n_constraints;
    std::unique_ptr<HostDeviceBuffer<ConstraintBond>> constraint_bonds;
    std::set<std::pair<int, int>> constrained_pairs;

    bool enabled() const {
        return n_constraints > 0;
    }
    bool contains(int ai, int aj) const {
        return constrained_pairs.count(std::minmax(ai, aj));
    }
};

class ConstraintForce {
   public:
    virtual ~ConstraintForce() = default;

    // Build constraints from the current context, remove constrained
    // bonded terms from the active topology, then initialize backend storage.
    void init(Context& ctx, const ParseResult& parsed);

    // Enforce constraints after an integration step. Implementations should
    // correct ctx.coords using ctx.xcoords as the reference frame.
    virtual void apply(Context& ctx, HostDeviceBuffer<coord_t>& xcoords) = 0;

    // Enforce constraints for a fresh start before the first MD step.
    // Implementations should make coordinates and velocities consistent with
    // the constrained geometry.
    virtual void initial_constraint(Context& ctx) = 0;

    /*
    Now we support mixed constriants force(LINCS, SHAKE, SETTLE) for different cases.
    So the main ConstraintForce will seperate the bonds to different constraint forces.
    */
    virtual void init_from_bonds(Context& ctx, const std::vector<ConstraintBond>& bonds) = 0;

    // Release backend-specific resources. CPU implementations usually have
    // nothing to free; CUDA implementations should free owned device state.
    virtual void cleanup() {}

    // True when at least one active constraint was built.
    bool enabled() const { return data_.enabled(); }

    // Access the built constraint data for tests or backend integration.
    ConstraintData& data() { return data_; }
    const ConstraintData& data() const { return data_; }

   protected:
    ConstraintData data_;

    // Backend hook called by init() after common constraint construction.
    // Use this to upload or transform data for a concrete implementation.
    virtual void init_backend(Context& ctx) {}

   private:
    // Populate data_ from topology, molecule, and SHAKE settings in ctx.
    void build_constraints(Context& ctx, const ParseResult& parsed);
};
