#pragma once

#include <cmath>

#include "precision.h"

inline force_accum_t force_to_accum(real_t value) {
    return static_cast<force_accum_t>(llrint(value * k_force_scale));
}

inline real_t force_from_accum(force_accum_t value) {
    return static_cast<real_t>(static_cast<double>(value) * static_cast<double>(k_inv_force_scale));
}

inline void add_force(force_accum_t& target, real_t value) {
    target += force_to_accum(value);
}

inline energy_accum_t energy_to_accum(real_t value) {
    return static_cast<energy_accum_t>(llrint(value * k_energy_scale));
}

inline real_t energy_from_accum(energy_accum_t value) {
    return static_cast<real_t>(static_cast<double>(value) * static_cast<double>(k_inv_energy_scale));
}

inline void add_energy(energy_accum_t& target, real_t value) {
    target += energy_to_accum(value);
}