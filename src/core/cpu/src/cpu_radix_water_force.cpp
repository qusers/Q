#include "cpu_radix_water_force.h"

#include <math.h>

#include "constants.h"
#include "context.h"

void calc_radix_w_forces() {
    auto& ctx = Context::instance();
    auto &coords = ctx.coords->cpu_data_p;
    auto &dvelocities = ctx.dvelocities->cpu_data_p;
    auto *excluded = ctx.excluded->cpu_data_p;

    real_t b, db, ener, dv, fexp;
    coord_t dr;
    real_t shift;

    if (ctx.md.radial_force != 0) {
        shift = sqrt(Boltz * ctx.Tfree / ctx.md.radial_force);
    } else {
        shift = 0;
    }

    for (int i = ctx.n_atoms_solute; i < ctx.n_atoms; i += 3) {
        if (excluded[i]) continue;

        dr.x = coords[i].x - ctx.topo.solvent_center.x;
        dr.y = coords[i].y - ctx.topo.solvent_center.y;
        dr.z = coords[i].z - ctx.topo.solvent_center.z;
        b = sqrt(dr.x * dr.x + dr.y * dr.y + dr.z * dr.z);
        db = b - (ctx.topo.solvent_radius - shift);

        if (db > 0) {
            ener = 0.5 * ctx.md.radial_force * db * db - ctx.Dwmz;
            dv = ctx.md.radial_force * db / b;
        } else if (b > 0.0) {
            fexp = exp(ctx.awmz * db);
            ener = ctx.Dwmz * (fexp * fexp - 2 * fexp);
            dv = -2 * ctx.Dwmz * ctx.awmz * (fexp - fexp * fexp) / b;
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
