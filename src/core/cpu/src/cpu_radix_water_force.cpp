#include "cpu_radix_water_force.h"

#include <math.h>

#include "constants.h"
#include "context.h"

void calc_radix_w_forces() {
    auto& ctx = Context::instance();
    auto &coords = ctx.coords->cpu_data_p;
    auto &dvelocities = ctx.dvelocities->cpu_data_p;

    real_t b, db, ener, dv, fexp;
    coord_t dr;
    real_t shift;

    if (ctx.md.radial_force != 0) {
        shift = sqrt(Boltz * ctx.Tfree / ctx.md.radial_force);
    } else {
        shift = 0;
    }

    for (int i = ctx.n_atoms_solute; i < ctx.n_atoms; i += 3) {
        dr.x = coords[i].x - ctx.topo.solvent_center.x;
        dr.y = coords[i].y - ctx.topo.solvent_center.y;
        dr.z = coords[i].z - ctx.topo.solvent_center.z;
        b = sqrt(pow(dr.x, 2) + pow(dr.y, 2) + pow(dr.z, 2));
        db = b - (ctx.topo.solvent_radius - shift);

        if (db > 0) {
            ener = 0.5 * ctx.md.radial_force * pow(db, 2) - ctx.Dwmz;
            dv = ctx.md.radial_force * db / b;
        } else if (b > 0.0) {
            fexp = exp(ctx.awmz * db);
            ener = ctx.Dwmz * (pow(fexp, 2) - 2 * fexp);
            dv = -2 * ctx.Dwmz * ctx.awmz * (fexp - pow(fexp, 2)) / b;
        } else {
            dv = 0;
            ener = 0;
        }

        ctx.E_restraint.Uradx += ener;
        dvelocities[i].x += dv * dr.x;
        dvelocities[i].y += dv * dr.y;
        dvelocities[i].z += dv * dr.z;
    }
}
