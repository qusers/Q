#include "cpu_temperature.h"

#include "constants.h"
#include "context.h"
#include "math.h"

/* =============================================
 * == ENERGY & TEMPERATURE
 * =============================================
 */

void calc_temperature() {
    auto& ctx = Context::instance();
    printf("Ndegf = %f, Ndegfree = %f, n_excluded = %d, Ndegfree_solvent = %f, Ndegfree_solute = %f\n",
           ctx.Ndegf, ctx.Ndegfree, ctx.n_excluded, ctx.Ndegfree_solvent, ctx.Ndegfree_solute);
    ctx.Temp = 0;
    ctx.Tfree = 0;
    double Temp_solute = 0, Tfree_solute = 0, Texcl_solute = 0;
    double Tfree_solvent = 0, Temp_solvent = 0, Texcl_solvent = 0;
    double Ekinmax = 1000.0 * ctx.Ndegf * Boltz * ctx.md.temperature / 2.0 / ctx.n_atoms;
    double ener;
    double mass_i;

    ctx.Temp = 0;
    for (int i = 0; i < ctx.n_atoms_solute; i++) {
        mass_i = ctx.catypes[ctx.atypes[i].code - 1].m;
        ener = .5 * mass_i * (pow(ctx.velocities[i].x, 2) + pow(ctx.velocities[i].y, 2) + pow(ctx.velocities[i].z, 2));
        Temp_solute += ener;
        if (!ctx.excluded[i]) {
            Tfree_solute += ener;
        } else {
            Texcl_solute += ener;
        }
        if (ener > Ekinmax) {
            printf(">>> WARNING: hot atom %d: %f\n", i, ener / Boltz / 3);
        }
    }

    for (int i = ctx.n_atoms_solute; i < ctx.n_atoms; i++) {
        mass_i = ctx.catypes[ctx.atypes[i].code - 1].m;
        ener = .5 * mass_i * (pow(ctx.velocities[i].x, 2) + pow(ctx.velocities[i].y, 2) + pow(ctx.velocities[i].z, 2));
        Temp_solvent += ener;
        if (!ctx.excluded[i]) {
            Tfree_solvent += ener;
        } else {
            Texcl_solvent += ener;
        }
        if (ener > Ekinmax) {
            printf(">>> WARNING: hot atom %d: %f\n", i, ener / Boltz / 3);
        }
    }

    ctx.Tfree = Tfree_solute + Tfree_solvent;
    ctx.Temp = Temp_solute + Temp_solvent;

    ctx.E_total.Ukin = ctx.Temp;

    ctx.Temp = 2.0 * ctx.Temp / Boltz / ctx.Ndegf;
    ctx.Tfree = 2.0 * ctx.Tfree / Boltz / ctx.Ndegfree;

    if (ctx.separate_scaling) {
        if (Tfree_solvent != 0) ctx.Tscale_solvent = sqrt(1 + (ctx.dt / ctx.tau_T) * (ctx.md.temperature / Tfree_solvent - 1.0));
        if (Tfree_solute != 0) ctx.Tscale_solute = sqrt(1 + (ctx.dt / ctx.tau_T) * (ctx.md.temperature / Tfree_solute - 1.0));
    } else {
        if (ctx.Tfree != 0) ctx.Tscale_solvent = sqrt(1 + (ctx.dt / ctx.tau_T) * (ctx.md.temperature / ctx.Tfree - 1.0));
        ctx.Tscale_solute = ctx.Tscale_solvent;
    }
    printf("Tscale = %f, tau_T = %f, Temp = %f, Tfree = %f\n", ctx.Tscale_solvent, ctx.tau_T, ctx.Temp, ctx.Tfree);
}