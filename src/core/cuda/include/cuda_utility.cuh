#pragma once

#include <math.h>

#include "constants.h"
#include "common/include/cuda_runtime_utility.h"
#include "common/include/precision.h"

__device__ inline real_t to_radians_device(real_t degrees) {
    return degrees * (static_cast<real_t>(q_fortran_pi) / 180.0);
}
