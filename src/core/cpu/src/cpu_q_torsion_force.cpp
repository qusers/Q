#include "cpu_q_torsion_force.h"

#include <math.h>

#include "context.h"
#include "cpu_force_accumulation.h"
#include "cpu_utils.h"

void calc_qtorsion_forces(int state) {
    auto& ctx = Context::instance();
    auto &coords = ctx.coords->cpu_data_p;
    auto &dvelocities = ctx.dvelocities->cpu_data_p;
    auto *lambdas = ctx.lambdas->cpu_data_p;
    int ic;
    int ai, aj, ak, al;
    coord_t rji, rjk, rkl, rnj, rnk, rki, rlj;
    coord_t di, dl, dpi, dpj, dpk, dpl;

    double bj2inv, bk2inv, bjinv, bkinv;
    double bj, bk, cos_phi, phi;
    double arg, dv, f1;
    double ener;
    energy_accum_t torsion = 0;

    for (int i = 0; i < ctx.n_qtorsions; i++) {
        ic = ctx.q_torsions[i + ctx.n_qtorsions * state].code;

        if (ic == 0) {
            continue;
        }

        ai = ctx.q_torsions[i + ctx.n_qtorsions * state].ai - 1;
        aj = ctx.q_torsions[i + ctx.n_qtorsions * state].aj - 1;
        ak = ctx.q_torsions[i + ctx.n_qtorsions * state].ak - 1;
        al = ctx.q_torsions[i + ctx.n_qtorsions * state].al - 1;

        rji.x = coords[ai].x - coords[aj].x;
        rji.y = coords[ai].y - coords[aj].y;
        rji.z = coords[ai].z - coords[aj].z;
        rjk.x = coords[ak].x - coords[aj].x;
        rjk.y = coords[ak].y - coords[aj].y;
        rjk.z = coords[ak].z - coords[aj].z;
        rkl.x = coords[al].x - coords[ak].x;
        rkl.y = coords[al].y - coords[ak].y;
        rkl.z = coords[al].z - coords[ak].z;
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
        dv = -ctx.q_ctorsions[ic].n * ctx.q_ctorsions[ic].k * sin(arg) * lambdas[state];

        // Forces
        f1 = sin(phi);
        if (fabs(f1) < k_singular_sin_epsilon) {
            f1 = copysign(k_singular_sin_epsilon, f1);
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
        add_energy(torsion, ener);

        add_force(dvelocities[ai].x, dv * dpi.x);
        add_force(dvelocities[ai].y, dv * dpi.y);
        add_force(dvelocities[ai].z, dv * dpi.z);

        add_force(dvelocities[aj].x, dv * dpj.x);
        add_force(dvelocities[aj].y, dv * dpj.y);
        add_force(dvelocities[aj].z, dv * dpj.z);

        add_force(dvelocities[ak].x, dv * dpk.x);
        add_force(dvelocities[ak].y, dv * dpk.y);
        add_force(dvelocities[ak].z, dv * dpk.z);

        add_force(dvelocities[al].x, dv * dpl.x);
        add_force(dvelocities[al].y, dv * dpl.y);
        add_force(dvelocities[al].z, dv * dpl.z);
    }

    ctx.energy.host()[ctx.energy.eq_index(ENERGY_FIXED_COUNT, state, EQ_BOND_TOR)] += torsion;
}
