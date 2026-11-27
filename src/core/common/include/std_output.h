#pragma once

#include "base_output.h"

class StdOutput : public BaseOutput {
   public:
    StdOutput() = default;
    ~StdOutput() override = default;
    OutputRequirements requirements(const Context& ctx, int iteration) const override;

   protected:
    void output_trajectory(Context& ctx, int iteration) override;
    void output_energy(Context& ctx, int iteration) override;
    void finish(Context& ctx) override;
};
