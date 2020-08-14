#ifndef __SOLVENT_H__
#define __SOLVENT_H__

void calc_nonbonded_ww_forces();
void calc_nonbonded_ww_forces_host();

void calc_nonbonded_pw_forces();
void calc_nonbonded_pw_forces_host();

struct calc_ww_t {
    dvel_t O;
    dvel_t H1;
    dvel_t H2;
};

struct calc_pw_t {
    dvel_t P;
    dvel_t W;
};

/* =============================================
 * == DEVICE
 * =============================================
 */

extern coord_t *P, *X;
extern dvel_t *DV_P, *DV;
extern calc_ww_t *WW_MAT;
extern calc_pw_t *PW_MAT;

// Constants pointers
extern ccharge_t *D_ccharges;
extern charge_t *D_charges;
extern catype_t *D_catypes;
extern atype_t *D_atypes;
extern p_atom_t *D_patoms;
extern bool *D_excluded;

// W-W interactions
__device__ void calc_ww_dvel_matrix_incr(int row, int column, double crg_ow, double crg_hw, double A_OO, double B_OO,
    coord_t *Xs, coord_t *Ys, double *Evdw, double *Ecoul, dvel_t *water_a, dvel_t *water_b);

__global__ void calc_ww_dvel_matrix(int n_waters, double crg_ow, double crg_hw, double A_OO, double B_OO,
    coord_t *X, double *Evdw, double *Ecoul, calc_ww_t *MAT);
    
__global__ void calc_ww_dvel_vector(int n_waters, dvel_t *DV, calc_ww_t *MAT);

// P-W interactions
__device__ void calc_pw_dvel_matrix_incr(int row, int pi, int column, int wj, int n_patoms,
    coord_t *Ps, coord_t *Xs, bool *excluded_s, double *Evdw, double *Ecoul, calc_pw_t *pw,
    ccharge_t *D_ccharges, charge_t *D_charges, catype_t *D_catypes, atype_t *D_atypes, p_atom_t *D_patoms);

__global__ void calc_pw_dvel_matrix(int n_patoms, int n_waters,
    coord_t *P, coord_t *X, double *Evdw, double *Ecoul, calc_pw_t *PW_MAT,
    ccharge_t *D_ccharges, charge_t *D_charges, catype_t *D_catypes, atype_t *D_atypes, p_atom_t *D_patoms, bool *D_excluded);

__global__ void calc_pw_dvel_vector_row(int n_patoms, int n_waters, dvel_t *DV_P, dvel_t *DV, calc_pw_t *PW_MAT);

__global__ void calc_pw_dvel_vector_column(int n_patoms, int n_waters, dvel_t *DV_P, dvel_t *DV, calc_pw_t *PW_MAT);


void clean_d_solvent();

#endif /* __SOLVENT_H__ */