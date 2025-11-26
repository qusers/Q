#ifndef CUDA_RESTRSEQ_FORCE_CUH
#define CUDA_RESTRSEQ_FORCE_CUH
#include "system.h"
__global__ void calc_restrseq_forces_kernel();

void calc_restrseq_forces_host();

void cleanup_restrseq_force();


#endif // CUDA_RESTRSEQ_FORCE_CUH