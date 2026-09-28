#pragma once

#if defined(QGPU_BACKEND_HIP)
#include <hip/hip_runtime.h>

using cudaError_t = hipError_t;

#define cudaSuccess           hipSuccess
#define cudaMalloc            hipMalloc
#define cudaFree              hipFree
#define cudaMemset            hipMemset
#define cudaMemcpy            hipMemcpy
#define cudaMemcpyHostToDevice hipMemcpyHostToDevice
#define cudaMemcpyDeviceToHost hipMemcpyDeviceToHost
#define cudaMemcpyDeviceToDevice hipMemcpyDeviceToDevice
#define cudaGetLastError      hipGetLastError
#define cudaGetErrorString    hipGetErrorString
#define cudaDeviceSynchronize hipDeviceSynchronize

#define cudaLaunchCooperativeKernel hipLaunchCooperativeKernel
#define cudaGetDevice hipGetDevice
#define cudaDeviceGetAttribute hipDeviceGetAttribute

#define cudaDevAttrMultiProcessorCount \
    hipDeviceAttributeMultiprocessorCount

#define cudaDevAttrCooperativeLaunch \
    hipDeviceAttributeCooperativeLaunch

#define cudaOccupancyMaxActiveBlocksPerMultiprocessor \
    hipOccupancyMaxActiveBlocksPerMultiprocessor

#else
#include <cuda_runtime.h>
#endif