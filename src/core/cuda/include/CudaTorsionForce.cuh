#ifndef CUDA_TORSION_FORCE_H
#define CUDA_TORSION_FORCE_H
#include "system.h"

__global__ void calc_torsion_forces_kernel(int start, int end, torsion_t* torsions, ctorsion_t* ctorsions, coord_t* coords, dvel_t* dvelocities, double* energy_sum);

double calc_torsion_forces_host(int start, int end);

void cleanup_torsion_force();

#endif