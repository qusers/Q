#include "cpu_q_improper_force.h"

#include <cmath>

#include "constants.h"
#include "context.h"
#include "cpu_utils.h"

void calc_qimproper_forces(int state) {
    auto& ctx = Context::instance();
    auto& coords = ctx.coords->cpu_data_p;
    auto& dvelocities = ctx.dvelocities->cpu_data_p;
    auto* lambdas = ctx.lambdas->cpu_data_p;

    for (int idx = 0; idx < ctx.n_qimpropers; idx++) {
        const q_improper_t& imp = ctx.q_impropers[idx + ctx.n_qimpropers * state];
        const int ic = imp.code;
        if (ic == 0) continue;

        const int aii = imp.ai - 1;
        const int aji = imp.aj - 1;
        const int aki = imp.ak - 1;
        const int ali = imp.al - 1;

        const coord_t ai = coords[aii];
        const coord_t aj = coords[aji];
        const coord_t ak = coords[aki];
        const coord_t al = coords[ali];

        coord_t rji = {ai.x - aj.x, ai.y - aj.y, ai.z - aj.z};
        coord_t rjk = {ak.x - aj.x, ak.y - aj.y, ak.z - aj.z};
        coord_t rkl = {al.x - ak.x, al.y - ak.y, al.z - ak.z};
        coord_t rnj = {
            rji.y * rjk.z - rji.z * rjk.y,
            rji.z * rjk.x - rji.x * rjk.z,
            rji.x * rjk.y - rji.y * rjk.x};
        coord_t rnk = {
            -rjk.y * rkl.z + rjk.z * rkl.y,
            -rjk.z * rkl.x + rjk.x * rkl.z,
            -rjk.x * rkl.y + rjk.y * rkl.x};

        const real_t bj2inv = static_cast<real_t>(1.0) /
                              (rnj.x * rnj.x + rnj.y * rnj.y + rnj.z * rnj.z);
        const real_t bk2inv = static_cast<real_t>(1.0) /
                              (rnk.x * rnk.x + rnk.y * rnk.y + rnk.z * rnk.z);
        const real_t bjinv = std::sqrt(bj2inv);
        const real_t bkinv = std::sqrt(bk2inv);

        real_t cos_phi = (rnj.x * rnk.x + rnj.y * rnk.y + rnj.z * rnk.z) *
                         (bjinv * bkinv);
        if (cos_phi > static_cast<real_t>(1.0)) cos_phi = static_cast<real_t>(1.0);
        if (cos_phi < static_cast<real_t>(-1.0)) cos_phi = static_cast<real_t>(-1.0);

        real_t phi = std::acos(cos_phi);
        const real_t sign =
            rjk.x * (rnj.y * rnk.z - rnj.z * rnk.y) +
            rjk.y * (rnj.z * rnk.x - rnj.x * rnk.z) +
            rjk.z * (rnj.x * rnk.y - rnj.y * rnk.x);
        if (sign < static_cast<real_t>(0.0)) phi = -phi;

        const q_cimproper_t& cimp = ctx.q_cimpropers[ic];
        real_t arg = phi - to_radians(cimp.phi0);
        arg -= static_cast<real_t>(2.0 * M_PI) *
               std::nearbyint(arg / static_cast<real_t>(2.0 * M_PI));
        real_t dv = cimp.k * arg;
        const real_t ener = static_cast<real_t>(0.5) * dv * arg;
        ctx.EQ_bond[state].Uimp += ener;
        dv *= lambdas[state];

        real_t f1 = std::sin(phi);
        if (std::fabs(f1) < k_singular_sin_epsilon) {
            f1 = k_singular_sin_epsilon;
        }
        f1 = static_cast<real_t>(-1.0) / f1;

        coord_t di = {
            f1 * (rnk.x * (bjinv * bkinv) - cos_phi * rnj.x * bj2inv),
            f1 * (rnk.y * (bjinv * bkinv) - cos_phi * rnj.y * bj2inv),
            f1 * (rnk.z * (bjinv * bkinv) - cos_phi * rnj.z * bj2inv)};
        coord_t dl = {
            f1 * (rnj.x * (bjinv * bkinv) - cos_phi * rnk.x * bk2inv),
            f1 * (rnj.y * (bjinv * bkinv) - cos_phi * rnk.y * bk2inv),
            f1 * (rnj.z * (bjinv * bkinv) - cos_phi * rnk.z * bk2inv)};

        coord_t rki = {rji.x - rjk.x, rji.y - rjk.y, rji.z - rjk.z};
        coord_t rlj = {-rjk.x - rkl.x, -rjk.y - rkl.y, -rjk.z - rkl.z};

        coord_t dpi = {
            rjk.y * di.z - rjk.z * di.y,
            rjk.z * di.x - rjk.x * di.z,
            rjk.x * di.y - rjk.y * di.x};
        coord_t dpj = {
            rki.y * di.z - rki.z * di.y + rkl.y * dl.z - rkl.z * dl.y,
            rki.z * di.x - rki.x * di.z + rkl.z * dl.x - rkl.x * dl.z,
            rki.x * di.y - rki.y * di.x + rkl.x * dl.y - rkl.y * dl.x};
        coord_t dpk = {
            rlj.y * dl.z - rlj.z * dl.y - rji.y * di.z + rji.z * di.y,
            rlj.z * dl.x - rlj.x * dl.z - rji.z * di.x + rji.x * di.z,
            rlj.x * dl.y - rlj.y * dl.x - rji.x * di.y + rji.y * di.x};
        coord_t dpl = {
            rjk.y * dl.z - rjk.z * dl.y,
            rjk.z * dl.x - rjk.x * dl.z,
            rjk.x * dl.y - rjk.y * dl.x};

        dvelocities[aii].x += dv * dpi.x;
        dvelocities[aii].y += dv * dpi.y;
        dvelocities[aii].z += dv * dpi.z;
        dvelocities[aji].x += dv * dpj.x;
        dvelocities[aji].y += dv * dpj.y;
        dvelocities[aji].z += dv * dpj.z;
        dvelocities[aki].x += dv * dpk.x;
        dvelocities[aki].y += dv * dpk.y;
        dvelocities[aki].z += dv * dpk.z;
        dvelocities[ali].x += dv * dpl.x;
        dvelocities[ali].y += dv * dpl.y;
        dvelocities[ali].z += dv * dpl.z;
    }
}
