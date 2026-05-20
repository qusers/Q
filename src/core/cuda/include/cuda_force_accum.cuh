#pragma once

#include "common/include/context.h"

__device__ __forceinline__ void atomic_add_force_component(force_accum_t* target, real_t value) {
    atomicAdd(target, static_cast<force_accum_t>(value));
}

__device__ __forceinline__ real_t force_accum_component_value(force_accum_t value) {
    return static_cast<real_t>(value);
}

__device__ __forceinline__ void atomic_add_force_component(force_fixed_storage_t* target, real_t value) {
    atomicAdd(target, force_fixed_to_storage(force_to_fixed(static_cast<double>(value))));
}

__device__ __forceinline__ real_t force_accum_component_value(force_fixed_storage_t value) {
    return fixed_to_force(force_fixed_from_storage(value));
}

inline cuda_dvel_t* cuda_force_accum_buffer(Context& ctx) {
    return ctx.fixed_dvelocities->gpu_data_p;
}
