#include "cpu_temperature.h"

#include "constants.h"
#include "context.h"
#include "cpu_force_accumulation.h"
#include "math.h"

void CpuTemperature::calc(Context& ctx) {
    const auto* masses = ctx.masses->cpu_data_p;
    auto& velocities = ctx.velocities->cpu_data_p;
    auto* excluded = ctx.excluded->cpu_data_p;

    energy_accum_t Temp_solute = 0, Tfree_solute = 0, Texcl_solute = 0;
    energy_accum_t Temp_solvent = 0, Tfree_solvent = 0, Texcl_solvent = 0;
    double Ekinmax = 1000.0 * data_.Ndegf * Boltz * ctx.md.temperature / 2.0 / ctx.n_atoms;  // Ekin_avg = Ndegf · kB · T / (2 · n_atoms), 1000 is tolerance factor
    for (int i = 0; i < ctx.n_atoms; i++) {
        double mass_i = masses[i];

        double v2 = velocities[i].x * velocities[i].x + velocities[i].y * velocities[i].y + velocities[i].z * velocities[i].z;

        double ener = 0.5 * mass_i * v2;

        if (ener > Ekinmax) {
            printf(">>> WARNING: hot atom %d: %f\n", i, ener / Boltz / 3);
        }

        if (i < ctx.n_atoms_solute) {
            // solute
            add_energy(Temp_solute, ener);
            if (!excluded[i]) {
                add_energy(Tfree_solute, ener);
            } else {
                add_energy(Texcl_solute, ener);
            }
        } else {
            // solvent
            add_energy(Temp_solvent, ener);
            if (!excluded[i]) {
                add_energy(Tfree_solvent, ener);
            } else {
                add_energy(Texcl_solvent, ener);
            }
        }
    }
    double Tfree_solute_value = energy_from_accum(Tfree_solute);
    double Tfree_solvent_value = energy_from_accum(Tfree_solvent);

    double* r = data_.results->cpu_data_p;

    double ukin_free = Tfree_solute_value + Tfree_solvent_value;
    double ukin = energy_from_accum(Temp_solute) + energy_from_accum(Temp_solvent);

    ctx.energy.data().Ukin = ukin;  // Now it is the system kinetic energy
    r[R_UKIN] = ukin;

    r[R_TEMP] = 2.0 * ukin / Boltz / data_.Ndegf;  // becomes a temperature Ekin = (Ndegf/2)·kB·T
    r[R_TFREE] = 2.0 * ukin_free / Boltz / data_.Ndegfree;

    if (data_.separate_scaling) {
        Tfree_solvent_value = 2.0 * Tfree_solvent_value / Boltz / data_.Ndegfree_solvent;
        Tfree_solute_value = 2.0 * Tfree_solute_value / Boltz / data_.Ndegfree_solute;

        //  Berendsen thermostat velocity-scaling factors
        if (Tfree_solvent_value != 0) {
            r[R_TSCALE_SLV] = sqrt(1 + (ctx.dt / data_.tau_T) * (ctx.md.temperature / Tfree_solvent_value - 1.0));
        }
        if (Tfree_solute_value != 0) {
            r[R_TSCALE_SOL] = sqrt(1 + (ctx.dt / data_.tau_T) * (ctx.md.temperature / Tfree_solute_value - 1.0));
        }
    } else {
        if (r[R_TFREE] != 0) {
            r[R_TSCALE_SLV] = sqrt(1 + (ctx.dt / data_.tau_T) * (ctx.md.temperature / r[R_TFREE] - 1.0));
        }
        r[R_TSCALE_SOL] = r[R_TSCALE_SLV];
    }

    printf("Tscale = %f, tau_T = %f, Temp = %f, Tfree = %f\n", r[R_TSCALE_SLV], data_.tau_T, r[R_TEMP], r[R_TFREE]);
}
