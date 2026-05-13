#pragma once

#include "common/include/precision.h"

void init_angle_force_kernel_data();
real_t calc_angle_forces_host(int start, int end);
void cleanup_angle_force();
