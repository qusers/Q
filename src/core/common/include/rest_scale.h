#pragma once
#include <cmath>

#include "cuda_runtime_utility.h"

// REST2 scale of the Q atoms: s = T0 / T_effective, 1 = off.
// Q-Q pairs and Q torsions are multiplied by s, Q-P and Q-W pairs by sqrt(s).
struct RestScale {
    double full = 1.0;  // s
    double half = 1.0;  // sqrt(s)

    // Scale of a pair that has n_q Q atoms (0, 1 or 2).
    HD double pair(int n_q) const { return n_q == 2 ? full : n_q == 1 ? half : 1.0; }
};

inline RestScale make_rest_scale(double s) { return {s, std::sqrt(s)}; }