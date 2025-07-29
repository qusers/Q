#ifndef __QATOMS_H__
#define __QATOMS_H__

void calc_nonbonded_qp_forces();
void calc_nonbonded_qp_forces_host();
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

struct calc_qp_t {
    dvel_t Q;
    dvel_t P;
};

/* =============================================
 * == DEVICE
 * =============================================
 */
extern coord_t *X, *W;
extern dvel_t *DV_X, *DV_W;
extern calc_qw_t *QW_MAT;
extern calc_qp_t *QP_MAT;
extern dvel_t *QQ_MAT;

// Constants pointers
extern ccharge_t *D_ccharges;
extern charge_t *D_charges;
extern catype_t *D_catypes;
extern atype_t *D_atypes;
extern p_atom_t *D_patoms;

extern q_catype_t *D_qcatypes;
extern q_atype_t *D_qatypes;
extern q_charge_t *D_qcharges;
extern q_atom_t *D_qatoms;
extern double *D_lambdas;
extern int *D_LJ_matrix;
extern bool *D_excluded;

// Q-W interactions
__device__ void calc_qw_dvel_matrix_incr(int row, int qi, int column, int n_lambdas, double crg_ow, double crg_hw, double A_O, double B_O,
    coord_t *Qs, coord_t *Ws, double *Evdw, double *Ecoul, calc_qw_t *qw,
    q_catype_t *D_qcatypes, q_atype_t *D_qatypes, q_charge_t *D_qcharges, q_atom_t *D_qatoms, double *D_lambdas, topo_t topo);

__global__ void calc_qw_dvel_matrix(int n_qatoms, int n_waters, int n_lambdas, double crg_ow, double crg_hw, double A_O, double B_O,
    coord_t *X, coord_t *W, double *Evdw, double *Ecoul, calc_qw_t *MAT,
    q_catype_t *D_qcatypes, q_atype_t *D_qatypes, q_charge_t *D_qcharges, q_atom_t *D_qatoms, double *D_lambdas, topo_t topo);

__global__ void calc_qw_dvel_vector_row(int n_qatoms, int n_waters, dvel_t *DV_X, dvel_t *DV_W, calc_qw_t *QW_MAT, q_atom_t *D_qatoms);

__global__ void calc_qw_dvel_vector_column(int n_qatoms, int n_waters, dvel_t *DV_X, dvel_t *DV_W, calc_qw_t *QW_MAT);

// Q-Q interactions

// Q-P interactions

__device__ void calc_qp_dvel_matrix_incr(int row, int qi, int column, int pj, int n_lambdas, int n_qatoms,
    coord_t *Qs, coord_t *Ps, int *LJs, bool *excluded_s, double *Evdw, double *Ecoul, calc_qp_t *qp,
    q_catype_t *D_qcatypes, q_atype_t *D_qatypes, q_charge_t *D_qcharges, p_atom_t *D_patoms, q_atom_t *D_qatoms, double *D_lambdas,
    catype_t *D_catypes, atype_t *D_atypes, ccharge_t *D_ccharges, charge_t *D_charges, topo_t topo);

__global__ void calc_qp_dvel_matrix(int n_qatoms, int n_patoms, int n_lambdas, int n_atoms_solute,
    coord_t *X, double *Evdw, double *Ecoul, calc_qp_t *QP_MAT,
    q_catype_t *D_qcatypes, q_atype_t *D_qatypes, q_charge_t *D_qcharges, p_atom_t *D_patoms, q_atom_t *D_qatoms, double *D_lambdas,
    int *D_LJ_matrix, bool *D_excluded, catype_t *D_catypes, atype_t *D_atypes, ccharge_t *D_ccharges, charge_t *D_charges, topo_t topo);

__global__ void calc_qp_dvel_vector_row(int n_qatoms, int n_patoms, dvel_t *DV_X, calc_qp_t *QP_MAT, q_atom_t *D_qatoms);

__global__ void calc_qp_dvel_vector_column(int n_qatoms, int n_patoms, dvel_t *DV_X, calc_qp_t *QP_MAT, p_atom_t *D_patoms);

__global__ void calc_energy_sum(int rows, int columns, double *Evdw_TOT, double *Ecoul_TOT, double *Evdw, double *Ecoul, bool upper_diagonal);

void clean_d_qatoms();

#endif /* __QATOMS_H__ */