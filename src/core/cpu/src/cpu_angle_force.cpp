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

    cangle_t cangle;
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

        const double rji_x = static_cast<double>(ai.x) - static_cast<double>(aj.x);
        const double rji_y = static_cast<double>(ai.y) - static_cast<double>(aj.y);
        const double rji_z = static_cast<double>(ai.z) - static_cast<double>(aj.z);

        const double rjk_x = static_cast<double>(ak.x) - static_cast<double>(aj.x);
        const double rjk_y = static_cast<double>(ak.y) - static_cast<double>(aj.y);
        const double rjk_z = static_cast<double>(ak.z) - static_cast<double>(aj.z);

        // Calculate inverse of norm of dist vector and their squares
        const double bji2inv = 1.0 / (rji_x * rji_x + rji_y * rji_y + rji_z * rji_z);
        const double bjk2inv = 1.0 / (rjk_x * rjk_x + rjk_y * rjk_y + rjk_z * rjk_z);
        const double bjiinv = sqrt(bji2inv);
        const double bjkinv = sqrt(bjk2inv);

        // Calculate cosine of angle and angle (th)
        double cos_th = (rji_x * rjk_x + rji_y * rjk_y + rji_z * rjk_z) * bjiinv * bjkinv;

        if (cos_th > 1.0) {
            cos_th = 1.0;
        } else if (cos_th < -1.0) {
            cos_th = -1.0;
        }

        const double th = acos(cos_th);

        const double dth = th - to_radians(cangle.th0);
        const double ener = 0.5 * cangle.kth * dth * dth;
        const double dv = cangle.kth * dth;

        double f1 = sin(th);
        const double sin_epsilon = static_cast<double>(k_singular_sin_epsilon);
        if (std::fabs(f1) < sin_epsilon) {
            // Avoid division by zero
            f1 = -1.0 / sin_epsilon;
        } else {
            f1 = -1.0 / f1;
        }

        // Update energies and forces

        add_energy(angle, ener);

        const double di_x = f1 * (rjk_x * bjiinv * bjkinv - cos_th * rji_x * bji2inv);
        const double di_y = f1 * (rjk_y * bjiinv * bjkinv - cos_th * rji_y * bji2inv);
        const double di_z = f1 * (rjk_z * bjiinv * bjkinv - cos_th * rji_z * bji2inv);

        const double dk_x = f1 * (rji_x * bjiinv * bjkinv - cos_th * rjk_x * bjk2inv);
        const double dk_y = f1 * (rji_y * bjiinv * bjkinv - cos_th * rjk_y * bjk2inv);
        const double dk_z = f1 * (rji_z * bjiinv * bjkinv - cos_th * rjk_z * bjk2inv);

        add_force(dvelocities[aii].x, dv * di_x);
        add_force(dvelocities[aii].y, dv * di_y);
        add_force(dvelocities[aii].z, dv * di_z);

        add_force(dvelocities[aki].x, dv * dk_x);
        add_force(dvelocities[aki].y, dv * dk_y);
        add_force(dvelocities[aki].z, dv * dk_z);

        add_force(dvelocities[aji].x, -dv * (di_x + dk_x));
        add_force(dvelocities[aji].y, -dv * (di_y + dk_y));
        add_force(dvelocities[aji].z, -dv * (di_z + dk_z));
    }

    return energy_from_accum(angle);
}
