#ifndef __SOLVENT_H__
#define __SOLVENT_H__

void calc_nonbonded_ww_forces();

__global__ void calc_ww_dvel_matrix(int n_waters, double crg_ow, double crg_hw, double A_OO, double B_OO,
    coord_t *d_X, dvel_t *d_DV, double *d_Evdw, double *d_Ecoul, dvel_t *d_MAT);
__global__ void calc_ww_dvel_vector(int n_waters, double crg_ow, double crg_hw, double A_OO, double B_OO,
    coord_t *d_X, dvel_t *d_DV, double *d_Evdw, double *d_Ecoul, dvel_t *d_MAT);

void calc_nonbonded_pw_forces();

#endif /* __SOLVENT_H__ */