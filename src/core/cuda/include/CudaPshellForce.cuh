#ifndef CUDA_PSHELL_FORCE_CUH
#define CUDA_PSHELL_FORCE_CUH
#include "system.h"
__global__ void calc_pshell_force_kernel(
    int n_atoms_solute,
    bool* shell,
    bool* excluded,
    coord_t* coords,
    coord_t* coords_top,
    double* ufix_energy,
    double* ushell_energy,
    dvel_t* dvelocities
);
void calc_pshell_forces_host();
void cleanup_pshell_force();

#endif  // CUDA_PSHELL_FORCE_CUH