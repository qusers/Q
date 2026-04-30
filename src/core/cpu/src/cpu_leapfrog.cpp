#include "cpu_leapfrog.h"

#include "context.h"
#include "cpu_shake.h"

void calc_leapfrog() {
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

        velocities[i].x = (velocities[i].x - dvelocities[i].x * ctx.dt * winv_i) * ctx.Tscale_solute;
        velocities[i].y = (velocities[i].y - dvelocities[i].y * ctx.dt * winv_i) * ctx.Tscale_solute;
        velocities[i].z = (velocities[i].z - dvelocities[i].z * ctx.dt * winv_i) * ctx.Tscale_solute;

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

        velocities[i].x = (velocities[i].x - dvelocities[i].x * ctx.dt * winv_i) * ctx.Tscale_solvent;
        velocities[i].y = (velocities[i].y - dvelocities[i].y * ctx.dt * winv_i) * ctx.Tscale_solvent;
        velocities[i].z = (velocities[i].z - dvelocities[i].z * ctx.dt * winv_i) * ctx.Tscale_solvent;

        xcoords[i].x = coords[i].x;
        xcoords[i].y = coords[i].y;
        xcoords[i].z = coords[i].z;

        coords[i].x += velocities[i].x * ctx.dt;
        coords[i].y += velocities[i].y * ctx.dt;
        coords[i].z += velocities[i].z * ctx.dt;
    }

    if (ctx.n_shake_constraints > 0) {
        calc_shake_constraints(coords, xcoords);
        for (int i = 0; i < ctx.n_atoms; i++) {
            velocities[i].x = (coords[i].x - xcoords[i].x) / ctx.dt;
            velocities[i].y = (coords[i].y - xcoords[i].y) / ctx.dt;
            velocities[i].z = (coords[i].z - xcoords[i].z) / ctx.dt;
        }
    }
}
