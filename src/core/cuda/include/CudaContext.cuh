#pragma once

#include <cuda_runtime.h>

#include "system.h"
#include "utils.h"

// CudaContext explicitly owns device buffers shared by CUDA kernels.
// Add/remove members to match the data you want to keep on the GPU.
class CudaContext {
   public:
    coord_t* d_coords = nullptr;
    dvel_t* d_dvelocities = nullptr;
    vel_t* d_velocities = nullptr;

    angle_t* d_angles = nullptr;
    cangle_t* d_cangles = nullptr;

    static CudaContext& instance() {
        static CudaContext ctx;
        return ctx;
    }
    void init();

    template <typename T>
    void sync_array_to_device(T* dst, const T* src, int count);

    template <typename T>
    void sync_array_to_host(T* dst, const T* src, int count);

    void sync_all_to_device();
    void sync_all_to_host();

   private:
    CudaContext() = default;

    void free();

    ~CudaContext() { free(); }
    CudaContext(const CudaContext&) = delete;
    CudaContext& operator=(const CudaContext&) = delete;
};
template <typename T>
void CudaContext::sync_array_to_device(T* dst, const T* src, int count) {
    check_cuda(cudaMemcpy(dst, src, count * sizeof(T), cudaMemcpyHostToDevice));
}
template <typename T>
void CudaContext::sync_array_to_host(T* dst, const T* src, int count) {
    check_cuda(cudaMemcpy(dst, src, count * sizeof(T), cudaMemcpyDeviceToHost));
}