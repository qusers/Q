#include "cpu_radix_water_force.h"

#include <math.h>

#include "constants.h"
#include "context.h"
#include "cpu_force_accumulation.h"

void calc_radix_w_forces() {
    auto& ctx = Context::instance();
    auto &coords = ctx.coords->cpu_data_p;
    auto &dvelocities = ctx.dvelocities->cpu_data_p;
    auto *excluded = ctx.excluded->cpu_data_p;

    energy_accum_t uradx = 0;

    double shift;
    if (ctx.md.radial_force != 0.0) {
        shift = sqrt(Boltz * ctx.Tfree / ctx.md.radial_force);
    } else {
        shift = 0.0;
    }

    for (int i = ctx.n_atoms_solute; i < ctx.n_atoms; i += 3) {
        if (excluded[i]) continue;

        const double dx = coords[i].x - ctx.topo.solvent_center.x;
        const double dy = coords[i].y - ctx.topo.solvent_center.y;
        const double dz = coords[i].z - ctx.topo.solvent_center.z;
        const double b = sqrt(dx * dx + dy * dy + dz * dz);
        const double db = b - (ctx.topo.solvent_radius - shift);

        double ener;
        double dv;
        if (db > 0.0) {
            ener = 0.5 * ctx.md.radial_force * db * db - ctx.Dwmz;
            dv = ctx.md.radial_force * db / b;
        } else if (b > 0.0) {
            const double fexp = exp(ctx.awmz * db);
            ener = ctx.Dwmz * (fexp * fexp - 2.0 * fexp);
            dv = -2.0 * ctx.Dwmz * ctx.awmz * (fexp - fexp * fexp) / b;
        } else {
            dv = 0.0;
            ener = 0.0;
        }

        add_energy(uradx, ener);
        add_force(dvelocities[i].x, dv * dx);
        add_force(dvelocities[i].y, dv * dy);
        add_force(dvelocities[i].z, dv * dz);
    }

    ctx.E_restraint.Uradx += energy_from_accum(uradx);
}
