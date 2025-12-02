#pragma once
#include "system.h"
void init_nonbonded_qw_force_kernel_data();

void calc_nonbonded_qw_forces_host_v2();

void cleanup_nonbonded_qw_force();