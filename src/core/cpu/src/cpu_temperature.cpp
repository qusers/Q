#include "cpu_temperature.h"

#include "constants.h"
#include "context.h"
#include "cpu_force_accumulation.h"
#include "math.h"

/* =============================================
 * == ENERGY & TEMPERATURE
 * =============================================
 */

void calc_temperature() {
    auto& ctx = Context::instance();
    auto &atypes = ctx.atypes->cpu_data_p;
    auto &catypes = ctx.catypes->cpu_data_p;
    auto &velocities = ctx.velocities->cpu_data_p;
    auto *excluded = ctx.excluded->cpu_data_p;
    ctx.Temp = 0;
    ctx.Tfree = 0;
    energy_accum_t Temp_solute = 0, Tfree_solute = 0, Texcl_solute = 0;
    energy_accum_t Tfree_solvent = 0, Temp_solvent = 0, Texcl_solvent = 0;
    double Ekinmax = 1000.0 * ctx.Ndegf * Boltz * ctx.md.temperature / 2.0 / ctx.n_atoms;
    double ener;
    double mass_i;

    ctx.Temp = 0;
    for (int i = 0; i < ctx.n_atoms_solute; i++) {
        mass_i = catypes[atypes[i].code - 1].m;
        ener = .5 * mass_i * (pow(velocities[i].x, 2) + pow(velocities[i].y, 2) + pow(velocities[i].z, 2));
        add_energy(Temp_solute, ener);
        if (!excluded[i]) {
            add_energy(Tfree_solute, ener);
        } else {
            add_energy(Texcl_solute, ener);
        }
        if (ener > Ekinmax) {
            printf(">>> WARNING: hot atom %d: %f\n", i, ener / Boltz / 3);
        }
    }

    for (int i = ctx.n_atoms_solute; i < ctx.n_atoms; i++) {
        mass_i = catypes[atypes[i].code - 1].m;
        ener = .5 * mass_i * (pow(velocities[i].x, 2) + pow(velocities[i].y, 2) + pow(velocities[i].z, 2));
        add_energy(Temp_solvent, ener);
        if (!excluded[i]) {
            add_energy(Tfree_solvent, ener);
        } else {
            add_energy(Texcl_solvent, ener);
        }
        if (ener > Ekinmax) {
            printf(">>> WARNING: hot atom %d: %f\n", i, ener / Boltz / 3);
        }
    }

    double Tfree_solute_value = energy_from_accum(Tfree_solute);
    double Tfree_solvent_value = energy_from_accum(Tfree_solvent);
    ctx.Tfree = Tfree_solute_value + Tfree_solvent_value;
    ctx.Temp = energy_from_accum(Temp_solute) + energy_from_accum(Temp_solvent);

    ctx.E_total.Ukin = ctx.Temp;

    ctx.Temp = 2.0 * ctx.Temp / Boltz / ctx.Ndegf;
    
    ctx.Tfree = 2.0 * ctx.Tfree / Boltz / ctx.Ndegfree;

    if (ctx.separate_scaling) {
        Tfree_solvent_value = 2.0 * Tfree_solvent_value / Boltz / ctx.Ndegfree_solvent;
        Tfree_solute_value = 2.0 * Tfree_solute_value / Boltz / ctx.Ndegfree_solute;
        if (Tfree_solvent_value != 0) ctx.Tscale_solvent = sqrt(1 + (ctx.dt / ctx.tau_T) * (ctx.md.temperature / Tfree_solvent_value - 1.0));
        if (Tfree_solute_value != 0) ctx.Tscale_solute = sqrt(1 + (ctx.dt / ctx.tau_T) * (ctx.md.temperature / Tfree_solute_value - 1.0));
    } else {
        if (ctx.Tfree != 0) ctx.Tscale_solvent = sqrt(1 + (ctx.dt / ctx.tau_T) * (ctx.md.temperature / ctx.Tfree - 1.0));
        ctx.Tscale_solute = ctx.Tscale_solvent;
    }

    printf("Tscale = %f, tau_T = %f, Temp = %f, Tfree = %f\n", ctx.Tscale_solvent, ctx.tau_T, ctx.Temp, ctx.Tfree);
}
