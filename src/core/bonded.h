#ifndef __BONDED_H__
#define __BONDED_H__

float calc_angle_forces(int start, int end);
float calc_bond_forces(int start, int end);
float calc_torsion_forces(int start, int end);
float calc_improper2_forces(int start, int end);

#endif /* __BONDED_H__ */