#pragma once

#include <fstream>
#include <string>
#include <vector>

#include "base_output.h"
#include "parse.h"

class NativeOutput : public BaseOutput {
   public:
    explicit NativeOutput(NativeOutputConfig config);
    ~NativeOutput() override = default;

    void init(Context& ctx) override;
    void finish(Context& ctx) override;
    void shutdown() override;

   protected:
    void output_trajectory(Context& ctx, int iteration) override;
    void output_energy(Context& ctx, int iteration) override;
    void output_restart(Context& ctx, int iteration) override;

   private:
    NativeOutputConfig config_;

    std::ofstream trajectory_stream_;
    std::ofstream energy_stream_;
    std::vector<int> trajectory_atoms_indices_;
    bool output_ready_ = false;

    void validate_config() const;
    void validate_runtime_config(Context& ctx) const;
    void ensure_initialized(Context& ctx);
    void init_trajectory_mask(Context& ctx);
    void open_trajectory(Context& ctx);
    void open_energy(Context& ctx);
    void write_trajectory_frame(Context& ctx);
    void write_energy_frame(Context& ctx);
    void write_restart_file(Context& ctx) const;
};
