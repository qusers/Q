#pragma once
#include "system.h"

void init_bond_force_kernel_data();
double calc_bond_forces_host(int start, int end);
void cleanup_bond_force();
