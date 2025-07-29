#ifndef __RESTRAINTS_H__
#define __RESTRAINTS_H__

void calc_radix_w_forces();
void calc_polx_w_forces(int iteration);
void calc_pshell_forces();

void calc_restrseq_forces();
void calc_restrdis_forces();
void calc_restrpos_forces();
void calc_restrang_forces();
void calc_restrwall_forces();

#endif  /* __RESTRAINTS_H__ */