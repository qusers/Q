#ifndef __QATOMS_H__
#define __QATOMS_H__

void calc_nonbonded_qp_forces();
void calc_nonbonded_qw_forces();
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


#endif /* __QATOMS_H__ */