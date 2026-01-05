#pragma once
#include "system.h"

void init_polx_water_force_kernel_data();

void calc_polx_water_forces_host(int iteration);

void cleanup_polx_water_force();
