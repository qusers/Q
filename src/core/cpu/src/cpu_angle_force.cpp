#include "cpu_angle_force.h"

#include <math.h>

#include "context.h"
#include "cpu_utils.h"

double calc_angle_forces(int start, int end) {
    auto& ctx = Context::instance();
    int aii, aji, aki;

    coord_t ai, aj, ak;
    coord_t rji, rjk;
    coord_t di, dk;

    double bji2inv, bjk2inv, bjiinv, bjkinv;
    cangle_t cangle;
    double cos_th, th, dth, dv, f1;
    double ener;
    double angle = 0;

    for (int i = start; i < end; i++) {
        aii = ctx.angles[i].ai - 1;
        aji = ctx.angles[i].aj - 1;
        aki = ctx.angles[i].ak - 1;
        ai = ctx.coords[aii];
        aj = ctx.coords[aji];
        ak = ctx.coords[aki];

        cangle = ctx.cangles[ctx.angles[i].code - 1];

        rji.x = ai.x - aj.x;
        rji.y = ai.y - aj.y;
        rji.z = ai.z - aj.z;

        rjk.x = ak.x - aj.x;
        rjk.y = ak.y - aj.y;
        rjk.z = ak.z - aj.z;

        // Calculate inverse of norm of dist vector and their squares
        bji2inv = 1 / (rji.x * rji.x + rji.y * rji.y + rji.z * rji.z);
        bjk2inv = 1 / (rjk.x * rjk.x + rjk.y * rjk.y + rjk.z * rjk.z);
        bjiinv = sqrt(bji2inv);
        bjkinv = sqrt(bjk2inv);

        // Calculate cosine of angle and angle (th)
        cos_th = (rji.x * rjk.x + rji.y * rjk.y + rji.z * rjk.z) * bjiinv * bjkinv;

        if (cos_th > 1) {
            cos_th = 1;
        } else if (cos_th < -1) {
            cos_th = -1;
        }

        th = acos(cos_th);

        dth = th - to_radians(cangle.th0);
        ener = .5 * cangle.kth * pow(dth, 2);
        dv = cangle.kth * dth;

        f1 = sin(th);
        if (std::fabs(f1) < 1.0E-12) {
            // Avoid division by zero
            f1 = -1.0E12;
        } else {
            f1 = -1.0 / f1;
        }

        // Update energies and forces

        angle += ener;

        di.x = f1 * (rjk.x * bjiinv * bjkinv - cos_th * rji.x * bji2inv);
        di.y = f1 * (rjk.y * bjiinv * bjkinv - cos_th * rji.y * bji2inv);
        di.z = f1 * (rjk.z * bjiinv * bjkinv - cos_th * rji.z * bji2inv);

        dk.x = f1 * (rji.x * bjiinv * bjkinv - cos_th * rjk.x * bjk2inv);
        dk.y = f1 * (rji.y * bjiinv * bjkinv - cos_th * rjk.y * bjk2inv);
        dk.z = f1 * (rji.z * bjiinv * bjkinv - cos_th * rjk.z * bjk2inv);

        ctx.dvelocities[aii].x += dv * di.x;
        ctx.dvelocities[aii].y += dv * di.y;
        ctx.dvelocities[aii].z += dv * di.z;

        ctx.dvelocities[aki].x += dv * dk.x;
        ctx.dvelocities[aki].y += dv * dk.y;
        ctx.dvelocities[aki].z += dv * dk.z;

        ctx.dvelocities[aji].x -= dv * (di.x + dk.x);
        ctx.dvelocities[aji].y -= dv * (di.y + dk.y);
        ctx.dvelocities[aji].z -= dv * (di.z + dk.z);

        // printf("ANGLE ener = %f\n", ener);
    }

    return angle;
}
