#ifndef CUDA_SHAKE_CONSTRAINTS_CUH
#define CUDA_SHAKE_CONSTRAINTS_CUH
#include "system.h"

__global__ void calc_shake_constraints_kernel();

void cleanup_shake_constraints();



#endif