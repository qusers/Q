#pragma once
#include "system.h"

void init_nonbonded_force_kernel_data();

std::pair<double, double> calc_nonbonded_force_host(int nx, int ny, int* x_idx_list, int* y_idx_list, bool symmetric);

void cleanup_nonbonded_force();