#include "cpu_bond_force.h"
#include "cpu_force_accumulation.h"

#include <math.h>

#include "context.h"

real_t calc_bond_forces(int start, int end) {
    auto& ctx = Context::instance();
    auto &bonds = ctx.bonds->cpu_data_p;
    auto &cbonds = ctx.cbonds->cpu_data_p;
    auto &coords = ctx.coords->cpu_data_p;
    auto &dvelocities = ctx.dvelocities->cpu_data_p;
    int aii, aji;
    coord_t ai, aj;
    cbond_t cbond;
    energy_accum_t bond = 0;

    for (int i = start; i < end; i++) {
        aii = bonds[i].ai - 1;
        aji = bonds[i].aj - 1;
        ai = coords[aii];
        aj = coords[aji];

        cbond = cbonds[bonds[i].code - 1];

        const double dx = static_cast<double>(aj.x) - static_cast<double>(ai.x);
        const double dy = static_cast<double>(aj.y) - static_cast<double>(ai.y);
        const double dz = static_cast<double>(aj.z) - static_cast<double>(ai.z);
        const double r2 = dx * dx + dy * dy + dz * dz;
        const double r = sqrt(r2);

        // Calculate energy
        const double dr = r - cbond.b0;
        const double ener = 0.5 * cbond.kb * dr * dr;

        add_energy(bond, ener);

        // Update forces
        const double ampl = cbond.kb * dr / r;

        add_force(dvelocities[aji].x, ampl * dx);
        add_force(dvelocities[aji].y, ampl * dy);
        add_force(dvelocities[aji].z, ampl * dz);

        add_force(dvelocities[aii].x, -ampl * dx);
        add_force(dvelocities[aii].y, -ampl * dy);
        add_force(dvelocities[aii].z, -ampl * dz);
    }

    return energy_from_accum(bond);
}
