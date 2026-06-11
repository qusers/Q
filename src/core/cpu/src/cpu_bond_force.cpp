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
    coord_t ai, aj, dx;
    cbond_t cbond;
    real_t dx2, dx1, ddx, ener, ampl;
    energy_accum_t bond = 0;

    for (int i = start; i < end; i++) {
        aii = bonds[i].ai - 1;
        aji = bonds[i].aj - 1;
        ai = coords[aii];
        aj = coords[aji];

        cbond = cbonds[bonds[i].code - 1];

        // Calculate distance vector, norm of distance vector
        dx.x = aj.x - ai.x;
        dx.y = aj.y - ai.y;
        dx.z = aj.z - ai.z;
        dx2 = dx.x * dx.x + dx.y * dx.y + dx.z * dx.z;
        dx1 = sqrt(dx2);

        // Calculate energy
        ddx = dx1 - cbond.b0;
        ener = .5 * cbond.kb * ddx * ddx;

        add_energy(bond, ener);

        // Update forces
        ampl = cbond.kb * ddx / dx1;

        add_force(dvelocities[aji].x, ampl * dx.x);
        add_force(dvelocities[aji].y, ampl * dx.y);
        add_force(dvelocities[aji].z, ampl * dx.z);

        add_force(dvelocities[aii].x, -ampl * dx.x);
        add_force(dvelocities[aii].y, -ampl * dx.y);
        add_force(dvelocities[aii].z, -ampl * dx.z);
    }

    return energy_from_accum(bond);
}
