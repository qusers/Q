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

    void init(Context& ctx);

    virtual void apply(Context& ctx) = 0;
    virtual void initial_shake(Context& ctx) = 0;
    virtual void cleanup() {}

    bool enabled() const { return data_.enabled(); }
    ShakeData& data() { return data_; }

   protected:
    ShakeData data_;
    virtual void init_backend(Context& ctx) {}

   private:
    void build_constraints(Context& ctx);
    void exclude_constrained_bonded_terms(Context& ctx);
};
