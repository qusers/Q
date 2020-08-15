#ifndef __QATOMS_H__
#define __QATOMS_H__

void calc_nonbonded_qp_forces();
void calc_nonbonded_qw_forces();
void calc_nonbonded_qw_forces_host();
void calc_nonbonded_qq_forces();

void calc_qangle_forces(int state);
void calc_qbond_forces(int state);
void calc_qtorsion_forces(int state);

struct calc_qw_t {
    dvel_t Q;
    dvel_t O;
    dvel_t H1;
    dvel_t H2;
};

/* =============================================
 * == DEVICE
 * =============================================
 */
extern coord_t *X, *W;
extern dvel_t *DV_X, *DV_W;
extern calc_qw_t *QW_MAT;

// Constants pointers
extern q_catype_t *D_qcatypes;
extern q_atype_t *D_qatypes;
extern q_charge_t *D_qcharges;
extern q_atom_t *D_qatoms;
extern double *D_lambdas;

// Q-W interactions
__device__ void calc_qw_dvel_matrix_incr(int row, int qi, int column, int n_lambdas, double crg_ow, double crg_hw, double A_O, double B_O,
    coord_t *Qs, coord_t *Ws, double *Evdw, double *Ecoul, calc_qw_t *qw,
    q_catype_t *D_qcatypes, q_atype_t *D_qatypes, q_charge_t *D_qcharges, q_atom_t *D_qatoms, double *D_lambdas);

__global__ void calc_qw_dvel_matrix(int n_qatoms, int n_waters, int n_lambdas, double crg_ow, double crg_hw, double A_O, double B_O,
    coord_t *Q, coord_t *W, double *Evdw, double *Ecoul, calc_qw_t *MAT,
    q_catype_t *D_qcatypes, q_atype_t *D_qatypes, q_charge_t *D_qcharges, q_atom_t *D_qatoms, double *D_lambdas);

__global__ void calc_qw_dvel_vector_row(int n_qatoms, int n_waters, dvel_t *DV_Q, dvel_t *DV_W, calc_qw_t *MAT, q_atom_t *D_qatoms);

__global__ void calc_qw_dvel_vector_column(int n_qatoms, int n_waters, dvel_t *DV_Q, dvel_t *DV_W, calc_qw_t *MAT);

// Q-Q interactions

// Q-P interactions

void clean_d_qatoms();

#endif /* __QATOMS_H__ */