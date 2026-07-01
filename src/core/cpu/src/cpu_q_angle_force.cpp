#include "cpu_q_angle_force.h"
#include "cpu_force_accumulation.h"

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
    energy_accum_t angle = 0;

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

        const double rji_x = coords[ai].x - coords[aj].x;
        const double rji_y = coords[ai].y - coords[aj].y;
        const double rji_z = coords[ai].z - coords[aj].z;

        const double rjk_x = coords[ak].x - coords[aj].x;
        const double rjk_y = coords[ak].y - coords[aj].y;
        const double rjk_z = coords[ak].z - coords[aj].z;

        const double bji2 = rji_x * rji_x + rji_y * rji_y + rji_z * rji_z;
        const double bjk2 = rjk_x * rjk_x + rjk_y * rjk_y + rjk_z * rjk_z;
        const double bji = sqrt(bji2);
        const double bjk = sqrt(bjk2);
        double cos_th = rji_x * rjk_x + rji_y * rjk_y + rji_z * rjk_z;
        cos_th /= (bji * bjk);
        if (cos_th > 1) {
            cos_th = 1;
        }
        if (cos_th < -1) {
            cos_th = -1;
        }
        const double th = acos(cos_th);
        const double dth = th - to_radians(ctx.q_cangles[ic].th0);
        const double ener = 0.5 * ctx.q_cangles[ic].kth * dth * dth;
        add_energy(angle, ener);

        const double dv = ctx.q_cangles[ic].kth * dth * lambdas[state];
        double f1 = sin(th);
        const double sin_epsilon = k_singular_sin_epsilon;
        if (fabs(f1) < sin_epsilon) {
            f1 = sin_epsilon;
        }
        f1 = -1.0 / f1;

        const double di_x = f1 * (rjk_x / (bji * bjk) - cos_th * rji_x / bji2);
        const double di_y = f1 * (rjk_y / (bji * bjk) - cos_th * rji_y / bji2);
        const double di_z = f1 * (rjk_z / (bji * bjk) - cos_th * rji_z / bji2);
        const double dk_x = f1 * (rji_x / (bji * bjk) - cos_th * rjk_x / bjk2);
        const double dk_y = f1 * (rji_y / (bji * bjk) - cos_th * rjk_y / bjk2);
        const double dk_z = f1 * (rji_z / (bji * bjk) - cos_th * rjk_z / bjk2);

        add_force(dvelocities[ai].x, dv * di_x);
        add_force(dvelocities[ai].y, dv * di_y);
        add_force(dvelocities[ai].z, dv * di_z);


        add_force(dvelocities[ak].x, dv * dk_x);
        add_force(dvelocities[ak].y, dv * dk_y);
        add_force(dvelocities[ak].z, dv * dk_z);


        add_force(dvelocities[aj].x, -dv * (di_x + dk_x));
        add_force(dvelocities[aj].y, -dv * (di_y + dk_y));
        add_force(dvelocities[aj].z, -dv * (di_z + dk_z));
    }

    ctx.energy.host()[ctx.energy.eq_index(ENERGY_FIXED_COUNT, state, EQ_BOND_ANGLE)] += angle;
}
