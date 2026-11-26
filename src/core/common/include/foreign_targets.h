#pragma once
#include "cuda_runtime_utility.h"
#include "rest_scale.h"

// Foreign-work targets per stage ([foreign_lambdas]); fixed so that the energy buffer keeps its size.
constexpr int MAX_FOREIGN_LAMBDAS = 8;

// The target windows of the current stage. Passed by value to the pair kernels.
// n == 0 means "do not evaluate", which is every step that does not write energies.
struct foreign_targets_t {
    int n = 0;
    double weight[MAX_FOREIGN_LAMBDAS][2] = {};  // lambda' of state 1 and 2
    double sc[MAX_FOREIGN_LAMBDAS][2] = {};      // (1 - lambda')^(1/6)
    RestScale rest[MAX_FOREIGN_LAMBDAS];         // REST2 scale of each target

    // Weight of a soft pair in target k: lambda' of its FEP state times the target's REST scale
    // for a pair with n_q Q atoms. The kernels add shape_weight * [V(sc') - V(sc)] to foreign_shape[k].
    HD double shape_weight(int k, int state, int n_q) const { return weight[k][state] * rest[k].pair(n_q); }
};
