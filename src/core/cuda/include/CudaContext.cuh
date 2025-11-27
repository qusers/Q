#pragma once

#include <cuda_runtime.h>

#include "system.h"
#include "utils.h"

class CudaContext {
   public:
    /*
    Common data
    */
    coord_t* d_coords = nullptr;
    dvel_t* d_dvelocities = nullptr;
    vel_t* d_velocities = nullptr;

    /*
    Used in CudaAngleForce.cu
    */
    angle_t* d_angles = nullptr;
    cangle_t* d_cangles = nullptr;

    /*
    Used in CudaBondForce.cu
    */
    bond_t* d_bonds = nullptr;
    cbond_t* d_cbonds = nullptr;

    /*
    Used in CudaImproper2Force.cu
    */
    improper_t* d_impropers = nullptr;
    cimproper_t* d_cimpropers = nullptr;

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