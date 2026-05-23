#pragma once

#include <math.h>

#include "common/include/precision.h"

template <typename Real>
QDYN_HOST_DEVICE inline Real q_softcore_pair_value(
    Real alpha_i,
    Real alpha_j,
    bool use_max_potential) {
    const Real eps = static_cast<Real>(1.0e-6);
    if (alpha_i <= eps && alpha_j <= eps) {
        return static_cast<Real>(0.0);
    }
    if (use_max_potential) {
        if (alpha_i > eps && alpha_j > eps) {
            return alpha_i < alpha_j ? alpha_i : alpha_j;
        }
        return alpha_i > alpha_j ? alpha_i : alpha_j;
    }
    return alpha_i > alpha_j ? alpha_i : alpha_j;
}

template <typename Real>
QDYN_HOST_DEVICE inline Real q_softcore_lookup_value(
    Real alpha,
    Real pair_a,
    Real pair_b,
    bool use_max_potential) {
    if (alpha == static_cast<Real>(0.0)) {
        return static_cast<Real>(0.0);
    }
    if (!use_max_potential) {
        return alpha;
    }

    const Real discriminant =
        pair_b * pair_b + static_cast<Real>(4.0) * alpha * pair_a;
    return (-pair_b + sqrt(discriminant)) /
           (static_cast<Real>(2.0) * alpha);
}

template <typename Real>
QDYN_HOST_DEVICE inline void q_softcore_lj(
    Real pair_a,
    Real pair_b,
    Real r6inv,
    Real softcore_lookup,
    Real* v_a,
    Real* v_b,
    Real* force_lj_term) {
    if (softcore_lookup == static_cast<Real>(0.0)) {
        *v_a = pair_a * r6inv * r6inv;
        *v_b = pair_b * r6inv;
        *force_lj_term = static_cast<Real>(12.0) * *v_a -
                         static_cast<Real>(6.0) * *v_b;
        return;
    }

    const Real r6_hc = static_cast<Real>(1.0) / r6inv;
    const Real r6_sc = r6_hc + softcore_lookup;
    const Real r6_sc_inv = static_cast<Real>(1.0) / r6_sc;
    *v_a = pair_a * r6_sc_inv * r6_sc_inv;
    *v_b = pair_b * r6_sc_inv;
    *force_lj_term = (static_cast<Real>(12.0) * *v_a -
                      static_cast<Real>(6.0) * *v_b) *
                     (r6_hc * r6_sc_inv);
}
