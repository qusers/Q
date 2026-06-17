#include "cpu_leapfrog.h"

#include "context.h"
#include "shake.h"
#include "cpu_force_accumulation.h"

void calc_leapfrog(Shake &shake) {
    auto& ctx = Context::instance();
    auto &atypes = ctx.atypes->cpu_data_p;
    auto &catypes = ctx.catypes->cpu_data_p;
    auto &coords = ctx.coords->cpu_data_p;
    auto &velocities = ctx.velocities->cpu_data_p;
    auto &dvelocities = ctx.dvelocities->cpu_data_p;
    auto *xcoords = ctx.xcoords->cpu_data_p;
    double mass_i;
    double winv_i;

    for (int i = 0; i < ctx.n_atoms_solute; i++) {
        mass_i = catypes[atypes[i].code - 1].m;
        winv_i = 1 / mass_i;

        const double fx = force_from_accum(dvelocities[i].x);
        const double fy = force_from_accum(dvelocities[i].y);
        const double fz = force_from_accum(dvelocities[i].z);

        velocities[i].x = (velocities[i].x - fx * ctx.dt * winv_i) * ctx.Tscale_solute;
        velocities[i].y = (velocities[i].y - fy * ctx.dt * winv_i) * ctx.Tscale_solute;
        velocities[i].z = (velocities[i].z - fz * ctx.dt * winv_i) * ctx.Tscale_solute;

        xcoords[i].x = coords[i].x;
        xcoords[i].y = coords[i].y;
        xcoords[i].z = coords[i].z;

        coords[i].x += velocities[i].x * ctx.dt;
        coords[i].y += velocities[i].y * ctx.dt;
        coords[i].z += velocities[i].z * ctx.dt;
    }

    for (int i = ctx.n_atoms_solute; i < ctx.n_atoms; i++) {
        mass_i = catypes[atypes[i].code - 1].m;
        winv_i = 1 / mass_i;

        const double fx = force_from_accum(dvelocities[i].x);
        const double fy = force_from_accum(dvelocities[i].y);
        const double fz = force_from_accum(dvelocities[i].z);

        velocities[i].x = (velocities[i].x - fx * ctx.dt * winv_i) * ctx.Tscale_solvent;
        velocities[i].y = (velocities[i].y - fy * ctx.dt * winv_i) * ctx.Tscale_solvent;
        velocities[i].z = (velocities[i].z - fz * ctx.dt * winv_i) * ctx.Tscale_solvent;

        xcoords[i].x = coords[i].x;
        xcoords[i].y = coords[i].y;
        xcoords[i].z = coords[i].z;

        coords[i].x += velocities[i].x * ctx.dt;
        coords[i].y += velocities[i].y * ctx.dt;
        coords[i].z += velocities[i].z * ctx.dt;
    }

    // do shake
    if (shake.data().n_constraints > 0) {
        shake.apply(ctx);
        for (int i = 0; i < ctx.n_atoms; i++) {
            velocities[i].x = (coords[i].x - xcoords[i].x) / ctx.dt;
            velocities[i].y = (coords[i].y - xcoords[i].y) / ctx.dt;
            velocities[i].z = (coords[i].z - xcoords[i].z) / ctx.dt;
        }
    }
}
