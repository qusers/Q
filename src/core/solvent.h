#ifndef __SOLVENT_H__
#define __SOLVENT_H__

void calc_nonbonded_ww_forces();

void calc_nonbonded_pw_forces();

struct calc_ww_t {
    dvel_t O;
    dvel_t H1;
    dvel_t H2;
};

struct calc_pw_t {
    dvel_t P;
    dvel_t W;
};


#endif /* __SOLVENT_H__ */