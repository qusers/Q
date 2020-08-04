#ifndef __SOLVENT_H__
#define __SOLVENT_H__

void calc_nonbonded_ww_forces();
void calc_nonbonded_ww_forces_host();

// Device pointers
extern coord_t *X;
extern dvel_t *MAT, *DV;

__device__ void calc_ww_dvel_matrix_incr(int row, int column, double crg_ow, double crg_hw, double A_OO, double B_OO,
    coord_t *Xs, coord_t *Ys, double *Evdw, double *Ecoul, dvel_t *water_a, dvel_t *water_b);

__global__ void calc_ww_dvel_matrix(int n_waters, double crg_ow, double crg_hw, double A_OO, double B_OO,
    coord_t *X, double *Evdw, double *Ecoul, dvel_t *MAT);
__global__ void calc_ww_dvel_vector_rows(int n_waters, dvel_t *DV, dvel_t *MAT);
__global__ void calc_ww_dvel_vector_columns(int n_waters, dvel_t *DV, dvel_t *MAT);

void calc_nonbonded_pw_forces();

#endif /* __SOLVENT_H__ */