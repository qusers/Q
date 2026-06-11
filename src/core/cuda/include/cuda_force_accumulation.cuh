#pragma once

#include <cuda_runtime.h>
#include <type_traits>

#include "precision.h"

template <typename T>
__device__ __forceinline__ force_accum_t force_to_accum(T value) {
    static_assert(std::is_same_v<T, float> || std::is_same_v<T, double>,
                  "force_to_accum only accepts float or double");

    return static_cast<force_accum_t>(llrint(value * k_force_scale));
}

__host__ __device__ __forceinline__ real_t force_from_accum(force_accum_t value) {
    return static_cast<real_t>(
        static_cast<double>(value) * static_cast<double>(k_inv_force_scale));
}

template <typename T>
__device__ __forceinline__ void atomic_add_force(force_accum_t* target, T value) {
    const force_accum_t fixed = force_to_accum(value);
    atomicAdd(
        reinterpret_cast<unsigned long long int*>(target),
        static_cast<unsigned long long int>(fixed));
}

template <typename T>
__device__ __forceinline__ energy_accum_t energy_to_accum(T value) {
    static_assert(std::is_same_v<T, float> || std::is_same_v<T, double>,
                  "energy_to_accum only accepts float or double");
    return static_cast<energy_accum_t>(llrint(value * k_energy_scale));
}

__host__ __device__ __forceinline__ real_t energy_from_accum(energy_accum_t value) {
    return static_cast<real_t>(
        static_cast<double>(value) * static_cast<double>(k_inv_energy_scale));
}

template <typename T>
__device__ __forceinline__ void atomic_add_energy(energy_accum_t* target, T value) {
    const energy_accum_t fixed = energy_to_accum(value);
    atomicAdd(
        reinterpret_cast<unsigned long long int*>(target),
        static_cast<unsigned long long int>(fixed));
}