// Idea for device pointers: allocate in system.cu, but have a struct of pointers to them in each
// forces calc file

#ifndef __NONBONDED_H__
#define __NONBONDED_H__

void calc_nonbonded_pp_forces();
void calc_nonbonded_pp_forces_host();

/* =============================================
 * == DEVICE
 * =============================================
 */

extern coord_t *X;
extern dvel_t *DV_X;
extern dvel_t *PP_MAT;

// Constants pointers
extern ccharge_t *D_ccharges;
extern charge_t *D_charges;
extern catype_t *D_catypes;
extern atype_t *D_atypes;
extern p_atom_t *D_patoms;
extern int *D_LJ_matrix;
extern bool *D_excluded;

// P-P interactions
__device__ void calc_pp_dvel_matrix_incr(int row, int pi, int column, int pj,
    coord_t *Xs, coord_t *Ys, int *LJs, bool *excluded_s, double *Evdw, double *Ecoul, dvel_t *patom_a, dvel_t *patom_b,
    ccharge_t *D_ccharges, charge_t *D_charges, catype_t *D_catypes, atype_t *D_atypes, p_atom_t *D_patoms, topo_t topo);

__global__ void calc_pp_dvel_matrix(int n_patoms, int n_atoms_solute,
    coord_t *X, double *Evdw, double *Ecoul, dvel_t *PP_MAT,
    ccharge_t *D_ccharges, charge_t *D_charges, catype_t *D_catypes, atype_t *D_atypes, p_atom_t *D_patoms, int *D_LJ_matrix, bool *D_excluded, topo_t topo);

__global__ void calc_pp_dvel_vector(int n_patoms, dvel_t *DV_X, dvel_t *PP_MAT, p_atom_t *D_patoms);


void clean_d_patoms();

__global__ void calc_energy_sum(int rows, int columns, double *Evdw_TOT, double *Ecoul_TOT, double *Evdw, double *Ecoul, bool upper_diagonal);

#endif /* __NONBONDED_H__ */