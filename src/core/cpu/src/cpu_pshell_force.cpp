#include "cpu_pshell_force.h"

#include <math.h>

#include "constants.h"
#include "context.h"

void calc_pshell_forces() {
    auto& ctx = Context::instance();
    auto &coords = ctx.coords->cpu_data_p;
    auto &dvelocities = ctx.dvelocities->cpu_data_p;
    auto *excluded = ctx.excluded->cpu_data_p;
    auto *shell = ctx.shell->cpu_data_p;

    coord_t dr;
    real_t k, r2, ener;

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
                ctx.E_restraint.Ufix += ener;
            }
            if (shell[i]) {
                ctx.E_restraint.Ushell += ener;
            }

            dvelocities[i].x += k * dr.x;
            dvelocities[i].y += k * dr.y;
            dvelocities[i].z += k * dr.z;
        }
    }
}
