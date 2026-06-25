#pragma once

#ifdef QDYN_SPFP
using real_t = float;
#else
using real_t = double;

#endif

using force_accum_t = long long;
using energy_accum_t = long long;

constexpr long long k_force_scale_i = 1ll << 40;
constexpr long long k_energy_scale_i = 1ll << 30;

#ifdef QDYN_SPFP
constexpr float k_force_scale = static_cast<float>(k_force_scale_i);
constexpr float k_energy_scale = static_cast<float>(k_energy_scale_i);
constexpr float k_inv_force_scale = static_cast<float>(1.0) / k_force_scale;
constexpr float k_inv_energy_scale = static_cast<float>(1.0) / k_energy_scale;
#else
constexpr double k_force_scale = static_cast<double>(k_force_scale_i);
constexpr double k_energy_scale = static_cast<double>(k_energy_scale_i);
constexpr double k_inv_force_scale = static_cast<double>(1.0) / k_force_scale;
constexpr double k_inv_energy_scale = static_cast<double>(1.0) / k_energy_scale;
#endif

constexpr double k_singular_sin_epsilon = 1.0e-12;