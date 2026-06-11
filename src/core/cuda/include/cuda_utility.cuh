#pragma once

#include <math.h>

#include "common/include/cuda_runtime_utility.h"
#include "common/include/precision.h"

__device__ inline double to_radians_device(double degrees) {
    return degrees * (M_PI / 180.0);
}
