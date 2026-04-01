#include "cpu_q_torsion_force.h"

#include <math.h>

#include "context.h"
#include "cpu_utils.h"

void calc_qtorsion_forces(int state) {
    auto& ctx = Context::instance();
    int ic;
    int ai, aj, ak, al;
    coord_t rji, rjk, rkl, rnj, rnk, rki, rlj;
    coord_t di, dl, dpi, dpj, dpk, dpl;

    double bj2inv, bk2inv, bjinv, bkinv;
    double bj, bk, cos_phi, phi;
    double arg, dv, f1;
    double ener;

    for (int i = 0; i < ctx.n_qtorsions; i++) {
        ic = ctx.q_torsions[i + ctx.n_qtorsions * state].code;

        if (ic == 0) {
            continue;
        }

        ai = ctx.q_torsions[i + ctx.n_qtorsions * state].ai - 1;
        aj = ctx.q_torsions[i + ctx.n_qtorsions * state].aj - 1;
        ak = ctx.q_torsions[i + ctx.n_qtorsions * state].ak - 1;
        al = ctx.q_torsions[i + ctx.n_qtorsions * state].al - 1;

        rji.x = ctx.coords[ai].x - ctx.coords[aj].x;
        rji.y = ctx.coords[ai].y - ctx.coords[aj].y;
        rji.z = ctx.coords[ai].z - ctx.coords[aj].z;
        rjk.x = ctx.coords[ak].x - ctx.coords[aj].x;
        rjk.y = ctx.coords[ak].y - ctx.coords[aj].y;
        rjk.z = ctx.coords[ak].z - ctx.coords[aj].z;
        rkl.x = ctx.coords[al].x - ctx.coords[ak].x;
        rkl.y = ctx.coords[al].y - ctx.coords[ak].y;
        rkl.z = ctx.coords[al].z - ctx.coords[ak].z;
        rnj.x = rji.y * rjk.z - rji.z * rjk.y;
        rnj.y = rji.z * rjk.x - rji.x * rjk.z;
        rnj.z = rji.x * rjk.y - rji.y * rjk.x;
        rnk.x = -rjk.y * rkl.z + rjk.z * rkl.y;
        rnk.y = -rjk.z * rkl.x + rjk.x * rkl.z;
        rnk.z = -rjk.x * rkl.y + rjk.y * rkl.x;

        bj = sqrt(pow(rnj.x, 2) + pow(rnj.y, 2) + pow(rnj.z, 2));
        bk = sqrt(pow(rnk.x, 2) + pow(rnk.y, 2) + pow(rnk.z, 2));
        cos_phi = (rnj.x * rnk.x + rnj.y * rnk.y + rnj.z * rnk.z) / (bj * bk);
        if (cos_phi > 1) {
            cos_phi = 1;
        }
        if (cos_phi < -1) {
            cos_phi = -1;
        }
        phi = acos(cos_phi);
        if (rjk.x * (rnj.y * rnk.z - rnj.z * rnk.y) + rjk.y * (rnj.z * rnk.x - rnj.x * rnk.z) +
                rjk.z * (rnj.x * rnk.y - rnj.y * rnk.x) <
            0) {
            phi = -phi;
        }

        bj2inv = 1 / (pow(rnj.x, 2) + pow(rnj.y, 2) + pow(rnj.z, 2));
        bk2inv = 1 / (pow(rnk.x, 2) + pow(rnk.y, 2) + pow(rnk.z, 2));
        bjinv = sqrt(bj2inv);
        bkinv = sqrt(bk2inv);

        // Energy
        arg = ctx.q_ctorsions[ic].n * phi - to_radians(ctx.q_ctorsions[ic].d);
        ener = ctx.q_ctorsions[ic].k * (1 + cos(arg));
        dv = -ctx.q_ctorsions[ic].n * ctx.q_ctorsions[ic].k * sin(arg) * ctx.lambdas[state];

        // Forces
        f1 = sin(phi);
        if (abs(f1) < 1E-12) {
            f1 = 1E-12;
        }
        f1 = -1 / f1;

        di.x = f1 * (rnk.x * (bjinv * bkinv) - cos_phi * rnj.x * bj2inv);
        di.y = f1 * (rnk.y * (bjinv * bkinv) - cos_phi * rnj.y * bj2inv);
        di.z = f1 * (rnk.z * (bjinv * bkinv) - cos_phi * rnj.z * bj2inv);
        dl.x = f1 * (rnj.x * (bjinv * bkinv) - cos_phi * rnk.x * bk2inv);
        dl.y = f1 * (rnj.y * (bjinv * bkinv) - cos_phi * rnk.y * bk2inv);
        dl.z = f1 * (rnj.z * (bjinv * bkinv) - cos_phi * rnk.z * bk2inv);

        rki.x = rji.x - rjk.x;
        rki.y = rji.y - rjk.y;
        rki.z = rji.z - rjk.z;
        rlj.x = -rjk.x - rkl.x;
        rlj.y = -rjk.y - rkl.y;
        rlj.z = -rjk.z - rkl.z;

        dpi.x = rjk.y * di.z - rjk.z * di.y;
        dpi.y = rjk.z * di.x - rjk.x * di.z;
        dpi.z = rjk.x * di.y - rjk.y * di.x;
        dpj.x = rki.y * di.z - rki.z * di.y + rkl.y * dl.z - rkl.z * dl.y;
        dpj.y = rki.z * di.x - rki.x * di.z + rkl.z * dl.x - rkl.x * dl.z;
        dpj.z = rki.x * di.y - rki.y * di.x + rkl.x * dl.y - rkl.y * dl.x;
        dpk.x = rlj.y * dl.z - rlj.z * dl.y - rji.y * di.z + rji.z * di.y;
        dpk.y = rlj.z * dl.x - rlj.x * dl.z - rji.z * di.x + rji.x * di.z;
        dpk.z = rlj.x * dl.y - rlj.y * dl.x - rji.x * di.y + rji.y * di.x;
        dpl.x = rjk.y * dl.z - rjk.z * dl.y;
        dpl.y = rjk.z * dl.x - rjk.x * dl.z;
        dpl.z = rjk.x * dl.y - rjk.y * dl.x;

        // Update energy and forces
        ctx.EQ_bond[state].Utor += ener;

        ctx.dvelocities[ai].x += dv * dpi.x;
        ctx.dvelocities[ai].y += dv * dpi.y;
        ctx.dvelocities[ai].z += dv * dpi.z;

        ctx.dvelocities[aj].x += dv * dpj.x;
        ctx.dvelocities[aj].y += dv * dpj.y;
        ctx.dvelocities[aj].z += dv * dpj.z;

        ctx.dvelocities[ak].x += dv * dpk.x;
        ctx.dvelocities[ak].y += dv * dpk.y;
        ctx.dvelocities[ak].z += dv * dpk.z;

        ctx.dvelocities[al].x += dv * dpl.x;
        ctx.dvelocities[al].y += dv * dpl.y;
        ctx.dvelocities[al].z += dv * dpl.z;
    }
}
