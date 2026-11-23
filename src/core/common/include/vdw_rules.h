// Defines vdW combination rules (geometric and arithmetic) and provides
// inline functions for computing vdW energy components for both rules.
#pragma once

#include <math.h>
#include <precision.h>
#include "constants.h"

#include "cuda_runtime_utility.h"

HD inline void calc_vdw_geometric(
    real_t ai_aii, real_t aj_aii, real_t ai_bii, real_t aj_bii,
    real_t r6, real_t* V_a, real_t* V_b) {
    *V_a = r6 * r6 * ai_aii * aj_aii;
    *V_b = r6 * ai_bii * aj_bii;
}

// Arithmetic rule: R*_ij = R*_i + R*_j, eps_ij = sqrt(eps_i * eps_j)
// Energy: V = sqrt(eps_ij) * R6^2 * r^-12 - 2*sqrt(eps_ij) * R6 * r^-6
// where R6 = (R*_i + R*_j)^6
// Parameters: For arithmetic rule, the same storage is reused:
//             ai_aii, aj_aii store R*_i, R*_j (vdW radius)
//             ai_bii, aj_bii store sqrt(eps_i), sqrt(eps_j) (after preprocessing)
//             r6 is 1/r^6
HD inline void calc_vdw_arithmetic(
    real_t Rstar_i, real_t Rstar_j, real_t sqrt_eps_i, real_t sqrt_eps_j,
    real_t r6, real_t* V_a, real_t* V_b) {
    real_t Rstar_ij = Rstar_i + Rstar_j;           // Arithmetic combination
    real_t sqrt_eps_ij = sqrt_eps_i * sqrt_eps_j;  // Geometric combination (already sqrt)

    // Compute R6 = (R*_ij)^6
    real_t R2 = Rstar_ij * Rstar_ij;
    real_t R6 = R2 * R2 * R2;

    *V_a = sqrt_eps_ij * R6 * R6 * r6 * r6;  // sqrt(eps_i * eps_j) * R^12 * r^-12
    *V_b = 2.0f * sqrt_eps_ij * R6 * r6;     // 2 * sqrt(eps_i * eps_j) * R^6 * r^-6
}

// LJ sigma of a pair from the per-atom parameters.
//   arithmetic: a = R*, so sigma = (R*_i + R*_j) / 2^(1/6)
//   geometric:  a = sqrt(C12), b = sqrt(C6), so sigma = (C12 / C6)^(1/6);
//               3 A when only one of C12, C6 is non-zero (GROMACS sc-sigma default).
HD inline real_t softcore_pair_sigma(int vdw_rule, real_t aii_i, real_t bii_i, real_t aii_j, real_t bii_j) {
    if (vdw_rule == VDW_GEOMETRIC) {
        const real_t c12 = aii_i * aii_j;
        const real_t c6 = bii_i * bii_j;
        if (c12 > 0 && c6 > 0) return pow(c12 / c6, static_cast<real_t>(1.0 / 6.0));
        if (c12 != 0 || c6 != 0) return static_cast<real_t>(3.0);
        return 0;
    }
    return (aii_i + aii_j) * static_cast<real_t>(0.8908987181403393);
}