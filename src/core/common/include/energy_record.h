#pragma once

#include <vector>

#include "base_output.h"

struct EnergyFrame {
    int step = 0;  // MD step; the final frame has step == md.steps
    double Upot = 0, Ukin = 0, Utot = 0;
    std::vector<double> q;  // per-state Q energy (EnergyData::eq_total), one per state
    std::vector<double> foreign_work;  // U(lambda') - U(lambda), one per foreign target
};

struct RunEnergies {
    std::vector<double> lambdas;      // λ of this run
    std::vector<double> foreign_lambdas;  // state-2 weight of each foreign target, as in StageSettings
    std::vector<EnergyFrame> frames;  // one per `energy` interval; empty when energy = 0
    EnergyFrame final;                // at the coordinates after the last step
};

class EnergyRecorder : public BaseOutput {
   public:
    explicit EnergyRecorder(RunEnergies& sink) : sink_(sink) {}
    void init(Context& ctx) override;
    void finish(Context& ctx) override;
    OutputRequirements requirements(const Context& ctx, int iteration) const override;

   protected:
    void output_trajectory(Context&, int) override {}
    void output_energy(Context& ctx, int iteration) override;

   private:
    RunEnergies& sink_;
    static EnergyFrame snapshot(Context& ctx, int step);
};