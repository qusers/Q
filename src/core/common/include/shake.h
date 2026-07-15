#pragma once
#include <memory>
#include <set>

#include "context.h"
#include "host_device_buffer.h"
#include "md_types.h"

struct ShakeBond {
    int ai;
    int aj;
    double dist2;
};

struct ShakeData {
    int n_constraints = 0;
    std::unique_ptr<HostDeviceBuffer<int>> mol_n_shakes;
    std::unique_ptr<HostDeviceBuffer<ShakeBond>> shake_bonds;
    std::set<std::pair<int, int>> constrained_pairs;

    bool enabled() const {
        return n_constraints > 0;
    }
    bool contains(int ai, int aj) const {
        return constrained_pairs.count(std::minmax(ai, aj));
    }
};

class Shake {
   public:
    virtual ~Shake() = default;

    // Build SHAKE constraints from the current context, remove constrained
    // bonded terms from the active topology, then initialize backend storage.
    void init(Context& ctx, const ParseResult& parsed);

    // Enforce constraints after an integration step. Implementations should
    // correct ctx.coords using ctx.xcoords as the reference frame.
    virtual void apply(Context &ctx, HostDeviceBuffer<coord_t>& xcoords) = 0;

    // Enforce constraints for a fresh start before the first MD step.
    // Implementations should make coordinates and velocities consistent with
    // the constrained geometry.
    virtual void initial_shake(Context& ctx, HostDeviceBuffer<coord_t>& xcoords_buffer) = 0;

    // Release backend-specific resources. CPU implementations usually have
    // nothing to free; CUDA implementations should free owned device state.
    virtual void cleanup() {}

    // True when at least one active SHAKE constraint was built.
    bool enabled() const { return data_.enabled(); }

    // Access the built constraint data for tests or backend integration.
    ShakeData& data() { return data_; }
    const ShakeData& data() const { return data_; }

   protected:
    ShakeData data_;

    // Backend hook called by init() after common constraint construction.
    // Use this to upload or transform data for a concrete implementation.
    virtual void init_backend(Context& ctx) {}

   private:
    // Populate data_ from topology, molecule, and SHAKE settings in ctx.
    void build_constraints(Context& ctx, const ParseResult& parsed);

    // Remove bonded force terms whose distances are handled by SHAKE.
    void exclude_constrained_bonded_terms(Context& ctx);
};
