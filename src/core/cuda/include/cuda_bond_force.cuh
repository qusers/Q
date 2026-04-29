#pragma once

#include "common/include/precision.h"

void init_bond_force_kernel_data();
real_t calc_bond_forces_host(int start, int end);
void cleanup_bond_force();
