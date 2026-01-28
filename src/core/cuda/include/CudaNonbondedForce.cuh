#pragma once
#include "system.h"

void init_nonbonded_force_kernel_data();

std::pair<double, double> calc_nonbonded_force_host(
    int nx, 
    int ny, 
    int* x_idx_list, 
    int* y_idx_list, 
    bool symmetric, 

    const int* x_charges_types, 
    const int* y_charges_types,
    const ccharge_t* ccharges_table,

    const int* x_atypes_types, 
    const int* y_atypes_types,
    const catype_t* catypes_table,
    const bool disable_water_h_lj = false,
    const double lambda = 1.0
);

void cleanup_nonbonded_force();
