#pragma once

#ifdef QDYN_SPFP
using real_t = float;
using nonbond_work_t = float;
using force_accum_t = float;
#else
using real_t = double;
using nonbond_work_t = double;
using force_accum_t = double;
#endif

using energy_accum_t = double;
using constraint_work_t = double;

#ifdef QDYN_SPFP
constexpr double k_singular_sin_epsilon = 1.0e-6;
#else
constexpr double k_singular_sin_epsilon = 1.0e-12;
#endif
