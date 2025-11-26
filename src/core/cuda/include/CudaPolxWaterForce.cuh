#ifndef __CUDA_POLX_WATER_FORCE_CUH__
#define __CUDA_POLX_WATER_FORCE_CUH__

#include "system.h"

void calc_polx_water_forces_host(int iteration);

void cleanup_polx_water_force();

#endif