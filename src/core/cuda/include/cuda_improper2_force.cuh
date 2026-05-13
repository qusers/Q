#pragma once

#include "common/include/precision.h"

void init_improper2_force_kernel_data();
real_t calc_improper2_forces_host(int start, int end);
void cleanup_improper2_force();
