// Defines vdW combination rules (geometric and arithmetic) and provides
// inline functions for computing vdW energy components for both rules.
#ifndef __VDW_RULES_H__
#define __VDW_RULES_H__

#include <math.h>

// Handle CUDA keywords for non-CUDA compilation
#ifndef __CUDACC__
#define __device__
#define __host__
#endif


// Geometric rule: A_ij = sqrt(A_i) * sqrt(A_j), B_ij = sqrt(B_i) * sqrt(B_j)
// Energy: V = A_ij * r^-12 - B_ij * r^-6
// Parameters: ai_aii, aj_aii are sqrt(A_i), sqrt(A_j)
//             ai_bii, aj_bii are sqrt(B_i), sqrt(B_j)
//             r6 is 1/r^6
__device__ __host__ inline void calc_vdw_geometric(
    double ai_aii, double aj_aii, double ai_bii, double aj_bii,
    double r6, double *V_a, double *V_b)
{
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
__device__ __host__ inline void calc_vdw_arithmetic(
    double Rstar_i, double Rstar_j, double sqrt_eps_i, double sqrt_eps_j,
    double r6, double *V_a, double *V_b)
{
    double Rstar_ij = Rstar_i + Rstar_j;         // Arithmetic combination
    double sqrt_eps_ij = sqrt_eps_i * sqrt_eps_j; // Geometric combination (already sqrt)

    // Compute R6 = (R*_ij)^6
    double R2 = Rstar_ij * Rstar_ij;
    double R6 = R2 * R2 * R2;

    *V_a = sqrt_eps_ij * R6 * R6 * r6 * r6;  // sqrt(eps_i * eps_j) * R^12 * r^-12
    *V_b = 2.0 * sqrt_eps_ij * R6 * r6;       // 2 * sqrt(eps_i * eps_j) * R^6 * r^-6
}

#endif /* __VDW_RULES_H__ */
