#include "cpu_angle_force.h"
#include "cpu_force_accumulation.h"

#include <math.h>

#include "context.h"
#include "cpu_utils.h"

real_t calc_angle_forces(int start, int end) {
    auto& ctx = Context::instance();
    auto &coords = ctx.coords->cpu_data_p;
    auto &dvelocities = ctx.dvelocities->cpu_data_p;
    int aii, aji, aki;

    coord_t ai, aj, ak;
    coord_t rji, rjk;
    coord_t di, dk;

    real_t bji2inv, bjk2inv, bjiinv, bjkinv;
    cangle_t cangle;
    real_t cos_th, th, dth, dv, f1;
    real_t ener;
    energy_accum_t angle = 0;

    auto &angles = ctx.angles->cpu_data_p;
    auto &cangles = ctx.cangles->cpu_data_p;

    for (int i = start; i < end; i++) {
        aii = angles[i].ai - 1;
        aji = angles[i].aj - 1;
        aki = angles[i].ak - 1;
        ai = coords[aii];
        aj = coords[aji];
        ak = coords[aki];

        cangle = cangles[angles[i].code - 1];

        rji.x = ai.x - aj.x;
        rji.y = ai.y - aj.y;
        rji.z = ai.z - aj.z;

        rjk.x = ak.x - aj.x;
        rjk.y = ak.y - aj.y;
        rjk.z = ak.z - aj.z;

        // Calculate inverse of norm of dist vector and their squares
        bji2inv = 1.0 / (rji.x * rji.x + rji.y * rji.y + rji.z * rji.z);
        bjk2inv = 1.0 / (rjk.x * rjk.x + rjk.y * rjk.y + rjk.z * rjk.z);
        bjiinv = sqrt(bji2inv);
        bjkinv = sqrt(bjk2inv);

        // Calculate cosine of angle and angle (th)
        cos_th = (rji.x * rjk.x + rji.y * rjk.y + rji.z * rjk.z) * bjiinv * bjkinv;

        if (cos_th > 1.0) {
            cos_th = 1.0;
        } else if (cos_th < -1.0) {
            cos_th = -1.0;
        }

        th = acos(cos_th);

        dth = th - to_radians(cangle.th0);
        ener = .5 * cangle.kth * pow(dth, 2);
        dv = cangle.kth * dth;

        f1 = sin(th);
        if (std::fabs(f1) < k_singular_sin_epsilon) {
            // Avoid division by zero
            f1 = -1.0 / k_singular_sin_epsilon;
        } else {
            f1 = -1.0 / f1;
        }

        // Update energies and forces

        add_energy(angle, ener);

        di.x = f1 * (rjk.x * bjiinv * bjkinv - cos_th * rji.x * bji2inv);
        di.y = f1 * (rjk.y * bjiinv * bjkinv - cos_th * rji.y * bji2inv);
        di.z = f1 * (rjk.z * bjiinv * bjkinv - cos_th * rji.z * bji2inv);

        dk.x = f1 * (rji.x * bjiinv * bjkinv - cos_th * rjk.x * bjk2inv);
        dk.y = f1 * (rji.y * bjiinv * bjkinv - cos_th * rjk.y * bjk2inv);
        dk.z = f1 * (rji.z * bjiinv * bjkinv - cos_th * rjk.z * bjk2inv);

        add_force(dvelocities[aii].x, dv * di.x);
        add_force(dvelocities[aii].y, dv * di.y);
        add_force(dvelocities[aii].z, dv * di.z);

        add_force(dvelocities[aki].x, dv * dk.x);
        add_force(dvelocities[aki].y, dv * dk.y);
        add_force(dvelocities[aki].z, dv * dk.z);

        add_force(dvelocities[aji].x, -dv * (di.x + dk.x));
        add_force(dvelocities[aji].y, -dv * (di.y + dk.y));
        add_force(dvelocities[aji].z, -dv * (di.z + dk.z));
    }

    return energy_from_accum(angle);
}
