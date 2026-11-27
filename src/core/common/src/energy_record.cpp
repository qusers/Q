#include "energy_record.h"

EnergyFrame EnergyRecorder::snapshot(Context& ctx, int step) {
    const EnergyData& e = ctx.energy.data();
    return {step, e.Upot, e.Ukin, e.Utot, e.eq_total, e.foreign_work};
}

void EnergyRecorder::init(Context& ctx) {
    sink_ = RunEnergies{};
    const double* l = ctx.lambdas->cpu_data_p;
    sink_.lambdas.assign(l, l + ctx.n_lambdas());
    for (int k = 0; k < ctx.foreign_targets.n; k++) sink_.foreign_lambdas.push_back(ctx.foreign_targets.weight[k][1]);
}

void EnergyRecorder::output_energy(Context& ctx, int iteration) {
    if (requirements(ctx, iteration).energy) sink_.frames.push_back(snapshot(ctx, iteration));
}

void EnergyRecorder::finish(Context& ctx) {
    sink_.final = snapshot(ctx, ctx.md.steps);  // Handler::finish_outputs() updated the totals first
}

// Same steps as the .en file (NativeOutput), so no extra GPU downloads.
OutputRequirements EnergyRecorder::requirements(const Context& ctx, int iteration) const {
    return {.energy = iteration > 0 && ctx.md.energy > 0 && iteration % ctx.md.energy == 0,
            .trajectory = false,
            .restart = false};
}

