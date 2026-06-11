#include "cpu_pshell_force.h"

#include <math.h>

#include "constants.h"
#include "context.h"
#include "cpu_force_accumulation.h"

void calc_pshell_forces() {
    auto& ctx = Context::instance();
    auto &coords = ctx.coords->cpu_data_p;
    auto &dvelocities = ctx.dvelocities->cpu_data_p;
    auto *excluded = ctx.excluded->cpu_data_p;
    auto *shell = ctx.shell->cpu_data_p;

    energy_accum_t ufix = 0;
    energy_accum_t ushell = 0;

    for (int i = 0; i < ctx.n_atoms_solute; i++) {
        if (shell[i] || excluded[i]) {
            const double k = excluded[i] ? k_fix : k_pshell;

            const double dx = static_cast<double>(coords[i].x) - static_cast<double>(ctx.coords_init->cpu_data_p[i].x);
            const double dy = static_cast<double>(coords[i].y) - static_cast<double>(ctx.coords_init->cpu_data_p[i].y);
            const double dz = static_cast<double>(coords[i].z) - static_cast<double>(ctx.coords_init->cpu_data_p[i].z);
            const double r2 = dx * dx + dy * dy + dz * dz;
            const double ener = 0.5 * k * r2;

            if (excluded[i]) {
                add_energy(ufix, ener);
            }
            if (shell[i]) {
                add_energy(ushell, ener);
            }

            add_force(dvelocities[i].x, k * dx);
            add_force(dvelocities[i].y, k * dy);
            add_force(dvelocities[i].z, k * dz);

        }
    }

    ctx.E_restraint.Ufix += energy_from_accum(ufix);
    ctx.E_restraint.Ushell += energy_from_accum(ushell);
}
