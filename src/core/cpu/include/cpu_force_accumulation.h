#pragma once

#include <cmath>
#include <type_traits>

#include "precision.h"

template <typename T>
inline force_accum_t force_to_accum(T value) {
    static_assert(std::is_same_v<T, float> || std::is_same_v<T, double>,
                  "force_to_accum only accepts float or double");

    return static_cast<force_accum_t>(llrint(value * k_force_scale));
}

inline double force_from_accum(force_accum_t value) {
    return static_cast<double>(static_cast<double>(value) * static_cast<double>(k_inv_force_scale));
}

template <typename T>
inline void add_force(force_accum_t& target, T value) {
    target += force_to_accum(value);
}

template <typename T>
inline energy_accum_t energy_to_accum(T value) {
    static_assert(std::is_same_v<T, float> || std::is_same_v<T, double>,
                  "energy_to_accum only accepts float or double");
    return static_cast<energy_accum_t>(llrint(value * k_energy_scale));
}

inline double energy_from_accum(energy_accum_t value) {
    return static_cast<double>(static_cast<double>(value) * static_cast<double>(k_inv_energy_scale));
}

template <typename T>
inline void add_energy(energy_accum_t& target, T value) {
    target += energy_to_accum(value);
}
