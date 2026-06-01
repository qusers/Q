#pragma once

#include <math.h>

#include "common/include/cuda_runtime_utility.h"
#include "common/include/precision.h"

__device__ inline real_t to_radians_device(real_t degrees) {
    return degrees * (M_PI / 180.0);
}
