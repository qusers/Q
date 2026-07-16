#include "cpu_integrator.h"

#include "cpu_force_accumulation.h"

void CpuIntegrator::step(Context& ctx) {
    auto& atypes = ctx.atypes->cpu_data_p;
    auto& catypes = ctx.catypes->cpu_data_p;
    auto& coords = ctx.coords->cpu_data_p;
    auto& velocities = ctx.velocities->cpu_data_p;
    auto& dvelocities = ctx.dvelocities->cpu_data_p;
    auto* xcoords = data_.xcoords->cpu_data_p;

    const double* temperature_results = temperature_->data().results->cpu_data_p;

    for (int i = 0; i < ctx.n_atoms; i++) {
        const double mass_i = catypes[atypes[i].code - 1].m;
        const double winv_i = 1.0 / mass_i;

        const double scale = (i < ctx.n_atoms_solute) ? temperature_results[R_TSCALE_SOL] : temperature_results[R_TSCALE_SLV];

        const double fx = force_from_accum(dvelocities[i].x);
        const double fy = force_from_accum(dvelocities[i].y);
        const double fz = force_from_accum(dvelocities[i].z);

        velocities[i].x = (velocities[i].x - fx * ctx.dt * winv_i) * scale;
        velocities[i].y = (velocities[i].y - fy * ctx.dt * winv_i) * scale;
        velocities[i].z = (velocities[i].z - fz * ctx.dt * winv_i) * scale;

        xcoords[i].x = coords[i].x;
        xcoords[i].y = coords[i].y;
        xcoords[i].z = coords[i].z;

        coords[i].x += velocities[i].x * ctx.dt;
        coords[i].y += velocities[i].y * ctx.dt;
        coords[i].z += velocities[i].z * ctx.dt;
    }

    if (shake_->data().n_constraints > 0) {
        shake_->apply(ctx, *data_.xcoords);
        for (int i = 0; i < ctx.n_atoms; i++) {
            velocities[i].x = (coords[i].x - xcoords[i].x) / ctx.dt;
            velocities[i].y = (coords[i].y - xcoords[i].y) / ctx.dt;
            velocities[i].z = (coords[i].z - xcoords[i].z) / ctx.dt;
        }
    }
}
