#ifndef CUDA_RESTRDIS_FORCE_CUH
#define CUDA_RESTRDIS_FORCE_CUH
#include "system.h"
__global__ void calc_restrdis_forces_kernel();

void calc_restrdis_forces_host();

void cleanup_restrdis_force();

#endif