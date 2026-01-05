#pragma once
#include "system.h"

void init_angle_force_kernel_data();
double calc_angle_forces_host(int start, int end);
void cleanup_angle_force();
