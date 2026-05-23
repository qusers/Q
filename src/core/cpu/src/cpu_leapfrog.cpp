#include "cpu_leapfrog.h"

#include "context.h"
#include "cpu_shake.h"

void calc_leapfrog() {
    auto& ctx = Context::instance();
    auto &coords = ctx.coords->cpu_data_p;
    auto &velocities = ctx.velocities->cpu_data_p;
    auto &dvelocities = ctx.dvelocities->cpu_data_p;
    auto *xcoords = ctx.xcoords->cpu_data_p;
    auto *winv = ctx.winv->cpu_data_p;
    real_t winv_i;

    for (int i = 0; i < ctx.n_atoms_solute; i++) {
        winv_i = winv[i];

        velocities[i].x = (velocities[i].x - dvelocities[i].x * winv_i * ctx.dt) * ctx.Tscale_solute;
        velocities[i].y = (velocities[i].y - dvelocities[i].y * winv_i * ctx.dt) * ctx.Tscale_solute;
        velocities[i].z = (velocities[i].z - dvelocities[i].z * winv_i * ctx.dt) * ctx.Tscale_solute;

        xcoords[i].x = coords[i].x;
        xcoords[i].y = coords[i].y;
        xcoords[i].z = coords[i].z;

        coords[i].x += velocities[i].x * ctx.dt;
        coords[i].y += velocities[i].y * ctx.dt;
        coords[i].z += velocities[i].z * ctx.dt;
    }

    for (int i = ctx.n_atoms_solute; i < ctx.n_atoms; i++) {
        winv_i = winv[i];

        velocities[i].x = (velocities[i].x - dvelocities[i].x * winv_i * ctx.dt) * ctx.Tscale_solvent;
        velocities[i].y = (velocities[i].y - dvelocities[i].y * winv_i * ctx.dt) * ctx.Tscale_solvent;
        velocities[i].z = (velocities[i].z - dvelocities[i].z * winv_i * ctx.dt) * ctx.Tscale_solvent;

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
