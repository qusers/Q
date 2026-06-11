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

    coord_t dr;
    real_t k, r2, ener;
    energy_accum_t ufix = 0;
    energy_accum_t ushell = 0;

    for (int i = 0; i < ctx.n_atoms_solute; i++) {
        if (shell[i] || excluded[i]) {
            if (excluded[i]) {
                k = k_fix;
            } else {
                k = k_pshell;
            }

            dr.x = coords[i].x - ctx.coords_init->cpu_data_p[i].x;
            dr.y = coords[i].y - ctx.coords_init->cpu_data_p[i].y;
            dr.z = coords[i].z - ctx.coords_init->cpu_data_p[i].z;
            r2 = pow(dr.x, 2) + pow(dr.y, 2) + pow(dr.z, 2);
            ener = 0.5 * k * r2;

            if (excluded[i]) {
                add_energy(ufix, ener);
            }
            if (shell[i]) {
                add_energy(ushell, ener);
            }

            add_force(dvelocities[i].x, k * dr.x);
            add_force(dvelocities[i].y, k * dr.y);
            add_force(dvelocities[i].z, k * dr.z);

        }
    }

    ctx.E_restraint.Ufix += energy_from_accum(ufix);
    ctx.E_restraint.Ushell += energy_from_accum(ushell);
}
