#pragma once

#include "base_output.h"

class StdOutput : public BaseOutput {
   public:
    StdOutput() = default;
    ~StdOutput() override = default;
    void output_final(Context& ctx, int iteration) override;

   protected:
    void output_trajectory(Context& ctx, int iteration) override;
    void output_energy(Context& ctx, int iteration) override;
};
