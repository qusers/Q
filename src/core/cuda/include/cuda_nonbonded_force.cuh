#pragma once

#include <utility>

#include "common/include/precision.h"

void init_nonbonded_force_kernel_data();

std::pair<real_t, real_t> calc_nonbonded_force_host(
    int nx, 
    int ny, 
    int* x_idx_list, 
    int* y_idx_list, 
    bool symmetric, 

    const int* x_charges_types, 
    const int* y_charges_types,
    const int* x_atypes_types, 
    const int* y_atypes_types,
    const bool disable_water_h_lj = false,
    const real_t lambda = 1.0,
    const bool softcore_enabled = false,
    const real_t* q_softcore_values = nullptr,
    const int* atom_to_qi = nullptr,
    const bool softcore_use_max_potential = false,
    const bool use_q_charge_product = false
);

void cleanup_nonbonded_force();
