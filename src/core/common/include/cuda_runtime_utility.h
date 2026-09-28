#pragma once

#include <cstdio>
#include <cstdlib>
#include "gpu_runtime_compat.h"

#if defined(__CUDACC__) || defined(__HIPCC__)
#define HD __host__ __device__
#else
#define HD
#endif

inline void check_cuda(cudaError_t status) {
    if (status != cudaSuccess) {
        std::printf(">>> FATAL: CUDA call failed with error code %d: %s\n", status, cudaGetErrorString(status));
        std::exit(EXIT_FAILURE);
    }
}

inline void check_cudaMalloc(void** dev_ptr, size_t size) {
    check_cuda(cudaMalloc(dev_ptr, size));
}
