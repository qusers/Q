#ifndef CUDA_TEMPERATURE_CUH
#define CUDA_TEMPERATURE_CUH
#include "system.h"
__global__ void calc_temperature_kernel(int n_atoms, int n_atoms_solute, atype_t* atypes, catype_t *catypes, vel_t *velocities, bool* excluded, double boltz, double ekinmax, 
double *Temp_solute, double *Tfree_solute, double *Texcl_solute, double *Temp_solvent, double *Tfree_solvent, double *Texcl_solvent);

void calc_temperature_host();

void cleanup_temperature();

#endif