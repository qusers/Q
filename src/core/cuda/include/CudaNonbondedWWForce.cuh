#pragma once
#include "system.h"
void init_nonbonded_ww_force_kernel_data();

void calc_nonbonded_ww_forces_host_v2();

void cleanup_nonbonded_ww_force();