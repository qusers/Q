#pragma once

#include <cuda_runtime.h>

#include <cstdio>
#include <cstdlib>

inline void check_cuda_impl(cudaError_t status, const char* file, int line, const char* expr) {
    if (status != cudaSuccess) {
        std::printf(">>> FATAL: CUDA call failed at %s:%d: %s -> error code %d: %s\n",
                    file, line, expr, status, cudaGetErrorString(status));
        std::exit(EXIT_FAILURE);
    }
}

#define check_cuda(status) check_cuda_impl((status), __FILE__, __LINE__, #status)

inline void check_cudaMalloc(void** dev_ptr, size_t size) {
    check_cuda(cudaMalloc(dev_ptr, size));
}
