#pragma once
#include "system.h"

void init_improper2_force_kernel_data();
double calc_improper2_forces_host(int start, int end);
void cleanup_improper2_force();
