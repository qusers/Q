#pragma once
#include <memory>

#include "host_device_buffer.h"
#include "md_types.h"
#include "context.h"


struct ShakeBond {
    int ai;
    int aj;
    real_t dist2;
};

struct ShakeData {
    int n_constraints = 0;
    std::unique_ptr<HostDeviceBuffer<int>> mol_n_shakes;
    std::unique_ptr<HostDeviceBuffer<ShakeBond>> shake_bonds;

    bool enabled() const {
        return n_constraints > 0;
    }
};

class Shake {
   public:
    virtual ~Shake() = default;

    // Build SHAKE constraints from the current context, remove constrained
    // bonded terms from the active topology, then initialize backend storage.
    void init(Context& ctx);

    // Enforce constraints after an integration step. Implementations should
    // correct ctx.coords using ctx.xcoords as the reference frame.
    virtual void apply(Context& ctx) = 0;

    // Enforce constraints for a fresh start before the first MD step.
    // Implementations should make coordinates and velocities consistent with
    // the constrained geometry.
    virtual void initial_shake(Context& ctx) = 0;

    // Release backend-specific resources. CPU implementations usually have
    // nothing to free; CUDA implementations should free owned device state.
    virtual void cleanup() {}

    // True when at least one active SHAKE constraint was built.
    bool enabled() const { return data_.enabled(); }

    // Access the built constraint data for tests or backend integration.
    ShakeData& data() { return data_; }

   protected:
    ShakeData data_;

    // Backend hook called by init() after common constraint construction.
    // Use this to upload or transform data for a concrete implementation.
    virtual void init_backend(Context& ctx) {}

   private:
    // Populate data_ from topology, molecule, and SHAKE settings in ctx.
    void build_constraints(Context& ctx);

    // Remove bonded force terms whose distances are handled by SHAKE.
    void exclude_constrained_bonded_terms(Context& ctx);
};
