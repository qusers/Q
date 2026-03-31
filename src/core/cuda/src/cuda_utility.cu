#include "cuda/include/cuda_utility.cuh"

#include <stdio.h>
#include <stdlib.h>

void check_cudaMalloc(void** devPtr, size_t size) {
    cudaError_t error = cudaMalloc((void**)devPtr, size);
    if (error != cudaSuccess) {
        printf(">>> FATAL: cudaMalloc failed with error code %d: %s\n", error, cudaGetErrorString(error));
        exit(EXIT_FAILURE);
    }
}

void check_cuda(cudaError_t status) {
    if (status != cudaSuccess) {
        printf(">>> FATAL: CUDA call failed with error code %d: %s\n", status, cudaGetErrorString(status));
        exit(EXIT_FAILURE);
    }
}
