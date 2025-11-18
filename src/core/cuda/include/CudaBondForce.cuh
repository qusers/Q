#ifndef CUDA_BOND_FORCE_H
#define CUDA_BOND_FORCE_H
#include "system.h"

__global__ void calc_bond_forces_kernel(int start, int end, bond_t* bonds, coord_t *coords, cbond_t *cbonds, dvel_t *dvelocities, double *energy_sum);
double calc_bond_forces_host(int start, int end);

#endif