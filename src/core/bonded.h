#ifndef __BONDED_H__
#define __BONDED_H__

double calc_angle_forces(int start, int end);
double calc_bond_forces(int start, int end);
double calc_torsion_forces(int start, int end);
double calc_improper2_forces(int start, int end);

#endif /* __BONDED_H__ */