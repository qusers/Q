#ifndef CUDA_RADIX_WATER_FORCE_CU
#define CUDA_RADIX_WATER_FORCE_CU
#include "system.h"


__global__ void calc_radix_water_forces_kernel(
    coord_t* coords, 
    double shift, 
    int n_atoms_solute,
    int n_atoms,
    topo_t topo,
    md_t md,
    double Dwmz,
    double awmz,
    dvel_t* dvelocities, 
    double* energy);
    
void calc_radix_water_forces_host();

void cleanup_radix_water_force();

#endif