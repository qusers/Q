#ifndef CUDA_SHAKE_CONSTRAINTS_CUH
#define CUDA_SHAKE_CONSTRAINTS_CUH
#include "system.h"

__global__ void calc_shake_constraints_kernel(
    int n_molecules,
    int *mol_n_shakes,
    shake_bond_t *shake_bonds,
    coord_t *coords,
    coord_t *xcoords,
    double *winv,
    int *total_iterations,
    int *mol_shake_offset
);

int calc_shake_constraints_host();

void cleanup_shake_constraints();



#endif