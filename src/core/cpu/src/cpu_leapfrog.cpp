#include "cpu_leapfrog.h"

#include "context.h"
#include "cpu_shake.h"

void calc_leapfrog() {
    auto& ctx = Context::instance();
    double mass_i;
    double winv_i;

    for (int i = 0; i < ctx.n_atoms_solute; i++) {
        mass_i = ctx.catypes[ctx.atypes[i].code - 1].m;
        winv_i = 1 / mass_i;

        ctx.velocities[i].x = (ctx.velocities[i].x - ctx.dvelocities[i].x * ctx.dt * winv_i) * ctx.Tscale_solute;
        ctx.velocities[i].y = (ctx.velocities[i].y - ctx.dvelocities[i].y * ctx.dt * winv_i) * ctx.Tscale_solute;
        ctx.velocities[i].z = (ctx.velocities[i].z - ctx.dvelocities[i].z * ctx.dt * winv_i) * ctx.Tscale_solute;

        ctx.xcoords[i].x = ctx.coords[i].x;
        ctx.xcoords[i].y = ctx.coords[i].y;
        ctx.xcoords[i].z = ctx.coords[i].z;

        ctx.coords[i].x += ctx.velocities[i].x * ctx.dt;
        ctx.coords[i].y += ctx.velocities[i].y * ctx.dt;
        ctx.coords[i].z += ctx.velocities[i].z * ctx.dt;
    }

    for (int i = ctx.n_atoms_solute; i < ctx.n_atoms; i++) {
        mass_i = ctx.catypes[ctx.atypes[i].code - 1].m;
        winv_i = 1 / mass_i;

        ctx.velocities[i].x = (ctx.velocities[i].x - ctx.dvelocities[i].x * ctx.dt * winv_i) * ctx.Tscale_solvent;
        ctx.velocities[i].y = (ctx.velocities[i].y - ctx.dvelocities[i].y * ctx.dt * winv_i) * ctx.Tscale_solvent;
        ctx.velocities[i].z = (ctx.velocities[i].z - ctx.dvelocities[i].z * ctx.dt * winv_i) * ctx.Tscale_solvent;

        ctx.xcoords[i].x = ctx.coords[i].x;
        ctx.xcoords[i].y = ctx.coords[i].y;
        ctx.xcoords[i].z = ctx.coords[i].z;

        ctx.coords[i].x += ctx.velocities[i].x * ctx.dt;
        ctx.coords[i].y += ctx.velocities[i].y * ctx.dt;
        ctx.coords[i].z += ctx.velocities[i].z * ctx.dt;
    }

    if (ctx.n_shake_constraints > 0) {
        calc_shake_constraints(ctx.coords.data(), ctx.xcoords.data());
        for (int i = 0; i < ctx.n_atoms; i++) {
            ctx.velocities[i].x = (ctx.coords[i].x - ctx.xcoords[i].x) / ctx.dt;
            ctx.velocities[i].y = (ctx.coords[i].y - ctx.xcoords[i].y) / ctx.dt;
            ctx.velocities[i].z = (ctx.coords[i].z - ctx.xcoords[i].z) / ctx.dt;
        }
    }
}
