#pragma once
#include <cmath>

#include "constants.h"
#include "cuda_runtime_utility.h"
#include "md_types.h"
#include "vdw_rules.h"
HD inline real_t2 calc_electrostatic(real_t qij, real_t coulomb_constant, real_t inv_dis) {
    if (qij == 0) return {0, 0};
    real_t vel = qij * coulomb_constant * inv_dis;  // k * qi * qj / r
    real_t dvel = -vel * inv_dis;                   /// -k * qi * qj / r2
    return {vel, dvel};
}

HD inline real_t2 calc_vdw(const real_t2& pair, real_t inv_dis) {
    if (pair.x == 0 && pair.y == 0) return {0, 0};
    real_t inv_dis3 = inv_dis * inv_dis * inv_dis;
    real_t inv_dis6 = inv_dis3 * inv_dis3;

    real_t V_a = pair.x * inv_dis6 * inv_dis6;  // precomputed A / r^12
    real_t V_b = pair.y * inv_dis6;             // precomputed B / r^6

    real_t vvdw = V_a - V_b;
    real_t dvvdw = (static_cast<real_t>(-12.0) * V_a + static_cast<real_t>(6.0) * V_b) * inv_dis;

    return {vvdw, dvvdw};
}

// Gapsys soft-core (Gapsys, Seeliger, de Groot, JCTC 2012), for r < rsc.
// The potential is a quadratic that matches the hard one's value, slope and curvature at rsc.
template <typename T>
HD inline Real2<T> calc_vdw_gapsys(const Real2<T>& pair, T dis, T rsc) {
    T inv = T(1.0) / rsc;
    T inv2 = inv * inv;
    T inv6 = inv2 * inv2 * inv2;

    T V_a = pair.x * inv6 * inv6;  // A / rsc^12
    T V_b = pair.y * inv6;         // B / rsc^6

    T a = (T(78.0) * V_a - T(21.0) * V_b) * inv2;
    T b = (T(168.0) * V_a - T(48.0) * V_b) * inv;
    T c = T(91.0) * V_a - T(28.0) * V_b;

    return {(a * dis - b) * dis + c, T(2.0) * a * dis - b};
}

template <typename T>
HD inline Real2<T> calc_electrostatic_gapsys(T qij, T coulomb_constant, T dis, T rsc) {
    T inv = T(1.0) / rsc;
    T qq = qij * coulomb_constant;

    T a = qq * inv * inv * inv;
    T b = T(3.0) * qq * inv * inv;
    T c = T(3.0) * qq * inv;

    return {(a * dis - b) * dis + c, T(2.0) * a * dis - b};
}

HD inline real_t2 combine_vdw(int vdw_rule, real_t aii_i, real_t bii_i, real_t aii_j, real_t bii_j) {
    real_t a, b;
    if (vdw_rule == VDW_GEOMETRIC) {
        calc_vdw_geometric(aii_i, aii_j, bii_i, bii_j, (real_t)1.0, &a, &b);
    } else {
        calc_vdw_arithmetic(aii_i, aii_j, bii_i, bii_j, (real_t)1.0, &a, &b);
    }
    return {a, b};
}

// Replaces el and vdw ({V, dV/dr}) of a soft pair where it is inside its soft-core radius.
HD inline void apply_gapsys(const softcore_pair_params_t& p, real_t sc, real_t sigma, real_t qij, real_t coulomb_constant, const real_t2& pair, real_t r, real_t2& el, real_t2& vdw) {
    const real_t rsc_lj = p.lj_scale * sigma * sc;
    if (r < rsc_lj && (pair.x != 0 || pair.y != 0)) {
        vdw = calc_vdw_gapsys(pair, r, rsc_lj);
    }

    const real_t rsc_q = p.alpha_q * (static_cast<real_t>(1.0) + p.sigma_q * fabs(qij)) * sc;
    if (r < rsc_q && qij != 0) {
        el = calc_electrostatic_gapsys(qij, coulomb_constant, r, rsc_q);
    }
}

// Energy of a pair at soft-core scale sc = (1 - lambda)^(1/6); sc <= 0 is the hard pair.
HD inline real_t gapsys_pair_energy(const softcore_pair_params_t& p, real_t sc, real_t sigma, real_t qij, real_t coulomb_constant, const real_t2& pair, real_t r) {
    const real_t inv_r = static_cast<real_t>(1.0) / r;
    real_t2 el = calc_electrostatic(qij, coulomb_constant, inv_r);
    real_t2 vdw = calc_vdw(pair, inv_r);
    if (sc > 0) apply_gapsys(p, sc, sigma, qij, coulomb_constant, pair, r, el, vdw);
    return el.x + vdw.x;
}

// V(r; sc_target) - V(r; sc_source) of a soft pair: what its energy changes by when only lambda changes.
// Zero outside both soft-core radii, where the two are the same hard pair.
HD inline real_t gapsys_shape_difference(const softcore_pair_params_t& p, real_t sc_source, real_t sc_target, real_t sigma, real_t qij,
                                         real_t coulomb_constant, const real_t2& pair, real_t r) {
    const real_t sc_max = sc_source > sc_target ? sc_source : sc_target;
    const real_t rsc_lj = p.lj_scale * sigma * sc_max;
    const real_t rsc_q = p.alpha_q * (static_cast<real_t>(1.0) + p.sigma_q * fabs(qij)) * sc_max;
    if (r <= 0 || (r >= rsc_lj && r >= rsc_q)) return 0;

    return gapsys_pair_energy(p, sc_target, sigma, qij, coulomb_constant, pair, r) -
           gapsys_pair_energy(p, sc_source, sigma, qij, coulomb_constant, pair, r);
}