#pragma once

#include <cuda_runtime.h>

#include <cstdio>
#include <cstdlib>

#include "qgpu_error.h"
#ifdef __CUDACC__
#define HD __host__ __device__

#else

#define HD
#endif

inline void check_cuda(cudaError_t status) {
    if (status != cudaSuccess) {
        throw CudaError(std::string("CUDA error ") + std::to_string(status) + ": " + cudaGetErrorString(status));
    }
}

// For destructors: never throw.
inline void check_cuda_noexcept(cudaError_t status) noexcept {
    if (status != cudaSuccess)
        std::fprintf(stderr, "CUDA error during cleanup: %s\n", cudaGetErrorString(status));
}

inline void check_cudaMalloc(void** dev_ptr, size_t size) {
    check_cuda(cudaMalloc(dev_ptr, size));
}
