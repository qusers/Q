#pragma once

#ifdef QDYN_SPFP
using real_t = float;
using nonbond_work_t = float;
#else
using real_t = double;
using nonbond_work_t = double;
#endif

using energy_accum_t = double;
using force_accum_t = double;
using constraint_work_t = double;
