#ifndef CUDA_ANGLE_FORCE_H
#define CUDA_ANGLE_FORCE_H
#include "system.h"
__global__ void calc_angle_forces_kernel(int start, int end, angle_t *angles, coord_t *coords, cangle_t *cangles, dvel_t *dvelocities, double *energy_sum);

double calc_angle_forces_host(int start, int end);
#endif