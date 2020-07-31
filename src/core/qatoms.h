#ifndef __QATOMS_H__
#define __QATOMS_H__

void calc_nonbonded_qp_qvdw_forces();
void calc_nonbonded_qp_forces();
void calc_nonbonded_qw_forces();
void calc_nonbonded_qq_forces();

void calc_qangle_forces(int state);
void calc_qbond_forces(int state);
void calc_qtorsion_forces(int state);

#endif /* __QATOMS_H__ */