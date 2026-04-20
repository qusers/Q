#include "cpu_pshell_force.h"

#include <math.h>

#include "constants.h"
#include "context.h"

void calc_pshell_forces() {
    auto& ctx = Context::instance();

    coord_t dr;
    double k, r2, ener;

    for (int i = 0; i < ctx.n_atoms_solute; i++) {
        if (ctx.shell[i] || ctx.excluded[i]) {
            if (ctx.excluded[i]) {
                k = k_fix;
            } else {
                k = k_pshell;
            }

            dr.x = ctx.coords[i].x - ctx.coords_init->cpu_data_p[i].x;
            dr.y = ctx.coords[i].y - ctx.coords_init->cpu_data_p[i].y;
            dr.z = ctx.coords[i].z - ctx.coords_init->cpu_data_p[i].z;
            r2 = pow(dr.x, 2) + pow(dr.y, 2) + pow(dr.z, 2);
            ener = 0.5 * k * r2;

            if (ctx.excluded[i]) {
                ctx.E_restraint.Ufix += ener;
            }
            if (ctx.shell[i]) {
                ctx.E_restraint.Ushell += ener;
            }

            ctx.dvelocities[i].x += k * dr.x;
            ctx.dvelocities[i].y += k * dr.y;
            ctx.dvelocities[i].z += k * dr.z;
        }
    }
}
