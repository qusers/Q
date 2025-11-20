#ifndef CUDA_UTILITY
#define CUDA_UTILITY

#include <math.h>

__device__ inline double to_radians_device(double degrees) {
    return degrees * (M_PI / 180.0);
}

#endif
