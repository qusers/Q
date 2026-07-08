#pragma once

void init_radix_water_force_kernel_data();
void calc_radix_water_forces_host(const double* temperature_results);
void cleanup_radix_water_force();
