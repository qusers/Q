#pragma once
#include "context.h"

struct OutputRequirements {
    bool energy = false;
    bool trajectory = false;
    bool restart = false;
};

class BaseOutput {
   public:
    explicit BaseOutput(int replica = 0) : replica_(replica) {}
    virtual ~BaseOutput() = default;
    virtual void init(Context& ctx) {}
    virtual void finish(Context& ctx) {}
    virtual void shutdown() {}
    void output(Context& ctx, int iteration) {
        output_trajectory(ctx, iteration);
        output_energy(ctx, iteration);
        output_restart(ctx, iteration);
    }
    virtual OutputRequirements requirements(
        const Context& ctx,
        int iteration) const = 0;

   protected:
    virtual void output_trajectory(Context& ctx, int iteration) = 0;
    virtual void output_energy(Context& ctx, int iteration) = 0;
    virtual void output_restart(Context& ctx, int iteration) {}
    int replica_ = 0;
};
