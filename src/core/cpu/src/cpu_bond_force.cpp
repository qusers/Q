#include "cpu_bond_force.h"

#include <math.h>

#include "context.h"

double calc_bond_forces(int start, int end) {
    auto& ctx = Context::instance();
    int aii, aji;
    coord_t ai, aj, dx;
    cbond_t cbond;
    double dx2, dx1, ddx, ener, ampl;
    double bond = 0;

    for (int i = start; i < end; i++) {
        aii = ctx.bonds[i].ai - 1;
        aji = ctx.bonds[i].aj - 1;
        ai = ctx.coords[aii];
        aj = ctx.coords[aji];

        cbond = ctx.cbonds[ctx.bonds[i].code - 1];

        // Calculate distance vector, norm of distance vector
        dx.x = aj.x - ai.x;
        dx.y = aj.y - ai.y;
        dx.z = aj.z - ai.z;
        dx2 = dx.x * dx.x + dx.y * dx.y + dx.z * dx.z;
        dx1 = sqrt(dx2);

        // Calculate energy
        ddx = dx1 - cbond.b0;
        ener = .5 * cbond.kb * ddx * ddx;

        bond += ener;

        // Update forces
        ampl = cbond.kb * ddx / dx1;

        ctx.dvelocities[aji].x += ampl * dx.x;
        ctx.dvelocities[aji].y += ampl * dx.y;
        ctx.dvelocities[aji].z += ampl * dx.z;

        ctx.dvelocities[aii].x -= ampl * dx.x;
        ctx.dvelocities[aii].y -= ampl * dx.y;
        ctx.dvelocities[aii].z -= ampl * dx.z;
    }

    return bond;
}
