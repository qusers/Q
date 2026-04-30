#include "cpu_q_angle_force.h"

#include <math.h>

#include "constants.h"
#include "context.h"
#include "cpu_utils.h"

void calc_qangle_forces(int state) {
    auto& ctx = Context::instance();
    auto &coords = ctx.coords->cpu_data_p;
    auto &dvelocities = ctx.dvelocities->cpu_data_p;
    auto *lambdas = ctx.lambdas->cpu_data_p;
    int ic;
    int ai, aj, ak;
    coord_t rji, rjk;
    double bji, bjk;
    double cos_th, th, dth, ener, dv, f1;
    coord_t di, dk;

    for (int i = 0; i < ctx.n_qangles; i++) {
        ic = ctx.q_angles[i + ctx.n_qangles * state].code - 1;

        // Skip if angle not present (code 0)
        // todo: Test it!!
        if (ic < 0) {
            continue;
        }

        ai = ctx.q_angles[i + ctx.n_qangles * state].ai - 1;
        aj = ctx.q_angles[i + ctx.n_qangles * state].aj - 1;
        ak = ctx.q_angles[i + ctx.n_qangles * state].ak - 1;

        rji.x = coords[ai].x - coords[aj].x;
        rji.y = coords[ai].y - coords[aj].y;
        rji.z = coords[ai].z - coords[aj].z;

        rjk.x = coords[ak].x - coords[aj].x;
        rjk.y = coords[ak].y - coords[aj].y;
        rjk.z = coords[ak].z - coords[aj].z;

        bji = sqrt(pow(rji.x, 2) + pow(rji.y, 2) + pow(rji.z, 2));
        bjk = sqrt(pow(rjk.x, 2) + pow(rjk.y, 2) + pow(rjk.z, 2));
        cos_th = rji.x * rjk.x + rji.y * rjk.y + rji.z * rjk.z;
        cos_th /= (bji * bjk);
        if (cos_th > 1) {
            cos_th = 1;
        }
        if (cos_th < -1) {
            cos_th = -1;
        }
        th = acos(cos_th);
        dth = th - to_radians(ctx.q_cangles[ic].th0);
        ener = .5 * ctx.q_cangles[ic].kth * pow(dth, 2);
        ctx.EQ_bond[state].Uangle += ener;

        dv = ctx.q_cangles[ic].kth * dth * lambdas[state];
        f1 = sin(th);
        if (abs(f1) < 1E-12) {
            f1 = 1E-12;
        }
        f1 = -1.0 / f1;

        di.x = f1 * (rjk.x / (bji * bjk) - cos_th * rji.x / pow(bji, 2));
        di.y = f1 * (rjk.y / (bji * bjk) - cos_th * rji.y / pow(bji, 2));
        di.z = f1 * (rjk.z / (bji * bjk) - cos_th * rji.z / pow(bji, 2));
        dk.x = f1 * (rji.x / (bji * bjk) - cos_th * rjk.x / pow(bjk, 2));
        dk.y = f1 * (rji.y / (bji * bjk) - cos_th * rjk.y / pow(bjk, 2));
        dk.z = f1 * (rji.z / (bji * bjk) - cos_th * rjk.z / pow(bjk, 2));

        dvelocities[ai].x += dv * di.x;
        dvelocities[ai].y += dv * di.y;
        dvelocities[ai].z += dv * di.z;
        dvelocities[ak].x += dv * dk.x;
        dvelocities[ak].y += dv * dk.y;
        dvelocities[ak].z += dv * dk.z;
        dvelocities[aj].x -= dv * (di.x + dk.x);
        dvelocities[aj].y -= dv * (di.y + dk.y);
        dvelocities[aj].z -= dv * (di.z + dk.z);
    }
}
