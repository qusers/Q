#pragma once
#include "system.h"

void init_nonbonded_qq_force_kernel_data();
void calc_nonbonded_qq_forces_host();

void cleanup_nonbonded_qq_force();
