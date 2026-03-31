#pragma once

#include <cuda_runtime.h>

#include <math.h>

__device__ inline double to_radians_device(double degrees) {
    return degrees * (M_PI / 180.0);
}

void check_cudaMalloc(void** devPtr, size_t size);
void check_cuda(cudaError_t status);
