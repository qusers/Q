#pragma once

#include <cstddef>
#include <cstring>
#include <stdexcept>

#include "cuda_runtime_utility.h"

template <typename T>
class HostDeviceBuffer {
   public:
    size_t length = 0;       // number of elements
    bool cpu_mem_flag = true;  // whether host memory exists
    bool gpu_mem_flag = true;  // whether device memory exists

    T* cpu_data_p = nullptr;  // CPU pointer
    T* gpu_data_p = nullptr;  // GPU pointer

    HostDeviceBuffer(int length, bool cpu_mem_flag = true, bool gpu_mem_flag = true);
    HostDeviceBuffer(unsigned int length, bool cpu_mem_flag = true, bool gpu_mem_flag = true);
    HostDeviceBuffer(size_t length, bool cpu_mem_flag = true, bool gpu_mem_flag = true);

    ~HostDeviceBuffer();

    void allocate();
    void deallocate();

    // copy CPU -> GPU
    void upload(T* buffer = nullptr);

    // copy GPU -> CPU
    void download(T* buffer = nullptr);
};

template <typename T>
HostDeviceBuffer<T>::HostDeviceBuffer(int length, bool cpu_mem_flag, bool gpu_mem_flag)
    : length(length), cpu_mem_flag(cpu_mem_flag), gpu_mem_flag(gpu_mem_flag) {
    allocate();
}

template <typename T>
HostDeviceBuffer<T>::HostDeviceBuffer(unsigned int length, bool cpu_mem_flag, bool gpu_mem_flag)
    : length(length), cpu_mem_flag(cpu_mem_flag), gpu_mem_flag(gpu_mem_flag) {
    allocate();
}

template <typename T>
HostDeviceBuffer<T>::HostDeviceBuffer(size_t length, bool cpu_mem_flag, bool gpu_mem_flag)
    : length(length), cpu_mem_flag(cpu_mem_flag), gpu_mem_flag(gpu_mem_flag) {
    allocate();
}

template <typename T>
HostDeviceBuffer<T>::~HostDeviceBuffer() {
    deallocate();
}

template <typename T>
void HostDeviceBuffer<T>::allocate() {
    if (length == 0) {
        return;
    }

    const size_t bytes = length * sizeof(T);
    if (cpu_mem_flag) {
        cpu_data_p = new T[length];
        std::memset(cpu_data_p, 0, bytes);
    }

    if (gpu_mem_flag) {
        check_cudaMalloc((void**)&gpu_data_p, bytes);
        check_cuda(cudaMemset(gpu_data_p, 0, bytes));
    }
}

template <typename T>
void HostDeviceBuffer<T>::deallocate() {
    if (cpu_data_p != nullptr) {
        delete[] cpu_data_p;
        cpu_data_p = nullptr;
    }

    if (gpu_data_p != nullptr) {
        check_cuda(cudaFree(gpu_data_p));
        gpu_data_p = nullptr;
    }

    length = 0;
}

//---------------------------------------------------------------------------------------------
// Upload: CPU -> GPU
//
// If buffer != nullptr, copy from buffer to device.
// Otherwise copy from internal cpu_data_p.
//
// Requires device memory.
//---------------------------------------------------------------------------------------------
template <typename T>
void HostDeviceBuffer<T>::upload(T* buffer) {
    if (length == 0) {
        return;
    }

    if (gpu_data_p == nullptr) {
        throw std::runtime_error("HostDeviceBuffer::upload failed: GPU memory is not allocated.");
    }

    T* src = buffer;
    if (src == nullptr) {
        src = cpu_data_p;
    }

    if (src == nullptr) {
        throw std::runtime_error("HostDeviceBuffer::upload failed: no CPU source buffer available.");
    }

    check_cuda(cudaMemcpy(gpu_data_p, src, length * sizeof(T), cudaMemcpyHostToDevice));
}

//---------------------------------------------------------------------------------------------
// Download: GPU -> CPU
//
// If buffer != nullptr, copy into buffer.
// Otherwise copy into internal cpu_data_p.
//
// Requires device memory.
//---------------------------------------------------------------------------------------------
template <typename T>
void HostDeviceBuffer<T>::download(T* buffer) {
    if (length == 0) {
        return;
    }

    if (gpu_data_p == nullptr) {
        throw std::runtime_error("HostDeviceBuffer::download failed: GPU memory is not allocated.");
    }

    T* dst = buffer;
    if (dst == nullptr) {
        dst = cpu_data_p;
    }

    if (dst == nullptr) {
        throw std::runtime_error("HostDeviceBuffer::download failed: no CPU destination buffer available.");
    }

    check_cuda(cudaMemcpy(dst, gpu_data_p, length * sizeof(T), cudaMemcpyDeviceToHost));
}
