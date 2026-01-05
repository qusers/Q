#pragma once
#include "system.h"

void init_torsion_force_kernel_data();
double calc_torsion_forces_host(int start, int end);

void cleanup_torsion_force();
