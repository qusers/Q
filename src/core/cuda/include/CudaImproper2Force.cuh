#ifndef CUDA_IMPROPER2_FORCE
#define CUDA_IMPROPER2_FORCE_H
#include "system.h"

__global__ void calc_improper2_forces_kernel(int start, int end, improper_t* impropers, cimproper_t* cimpropers, coord_t* coords, dvel_t* dvelocities, double* energy_sum);

double calc_improper2_forces_host(int start, int end);

void cleanup_improper2_force();

#endif