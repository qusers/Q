#pragma once

#include <limits>

#if defined(__CUDACC__)
#define QDYN_HOST_DEVICE __host__ __device__
#else
#define QDYN_HOST_DEVICE
#endif

#ifdef QDYN_SPFP
using real_t = float;
using nonbond_work_t = float;
using force_accum_t = float;
using energy_accum_t = float;
using constraint_work_t = double;
#else
using real_t = double;
using nonbond_work_t = double;
using force_accum_t = double;
using energy_accum_t = double;
using constraint_work_t = double;
#endif

using force_fixed_t = long long;
using force_fixed_storage_t = unsigned long long;

constexpr double k_force_fixed_scale = 1.0e8;
constexpr double k_force_fixed_max_scaled = 9.2233720368547758e18;

QDYN_HOST_DEVICE inline force_fixed_t force_to_fixed(double value) {
    const double scaled = value * k_force_fixed_scale;
#ifndef NDEBUG
    if (!(scaled <= k_force_fixed_max_scaled && scaled >= -k_force_fixed_max_scaled)) {
#if defined(__CUDA_ARCH__)
        asm("trap;");
#else
        __builtin_trap();
#endif
    }
#endif
    return scaled >= 0.0 ? static_cast<force_fixed_t>(scaled + 0.5) : static_cast<force_fixed_t>(scaled - 0.5);
}

QDYN_HOST_DEVICE inline force_fixed_storage_t force_fixed_to_storage(force_fixed_t value) {
    return static_cast<force_fixed_storage_t>(value);
}

QDYN_HOST_DEVICE inline force_fixed_t force_fixed_from_storage(force_fixed_storage_t value) {
    constexpr force_fixed_storage_t k_sign_bit = 1ULL << 63;
    if ((value & k_sign_bit) == 0) {
        return static_cast<force_fixed_t>(value);
    }
    if (value == k_sign_bit) {
        return (-9223372036854775807LL - 1LL);
    }
    return -static_cast<force_fixed_t>((~value) + 1ULL);
}

QDYN_HOST_DEVICE inline real_t fixed_to_force(force_fixed_t value) {
    return static_cast<real_t>(static_cast<double>(value) / k_force_fixed_scale);
}

#ifdef QDYN_SPFP
constexpr double k_singular_sin_epsilon = 1.0e-6;
#else
constexpr double k_singular_sin_epsilon = 1.0e-12;
#endif
