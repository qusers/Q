#pragma once

#include <string>
#include <vector>

#include "common/include/precision.h"
/* =============================================
 * == FROM MD FILE
 * =============================================
 */

struct md_t {
    // [MD]
    int steps;
    real_t stepsize;
    real_t temperature;
    char thermostat[40];
    real_t bath_coupling;
    int random_seed;
    real_t initial_temperature;
    bool shake_solvent;
    bool shake_solute;
    bool shake_hydrogens;
    bool lrf;
    bool charge_groups;
    // [cut-offs]
    real_t solute_solute;
    real_t solvent_solvent;
    real_t solute_solvent;
    real_t q_atom;
    // [sphere]
    real_t shell_radius;  // Note: this is for the pshell
    real_t shell_force;   // Note: this is for the pshell
    // [solvent]
    real_t radial_force;
    bool polarisation;
    real_t polarisation_force;
    // [intervals]
    int non_bond;
    int output;
    int energy;
    int trajectory;
    // [trajectory_atoms]
    // [lambdas]
    // [sequence_restraints]
    // [distance_restraints]
    // [angle_restraints]
    // [wall_restraints]
};

struct coord_t {
    real_t x;
    real_t y;
    real_t z;
};

struct bond_t {
    int ai;
    int aj;
    int code;
};

struct cbond_t {
    int code;
    real_t kb;
    real_t b0;
};

struct angle_t {
    int ai;
    int aj;
    int ak;
    int code;
};

struct cangle_t {
    int code;
    real_t kth;
    real_t th0;
};

struct torsion_t {
    int ai;
    int aj;
    int ak;
    int al;
    int code;
};

struct ctorsion_t {
    int code;
    real_t k;
    real_t n;
    real_t d;
    real_t paths;
};

struct improper_t {
    int ai;
    int aj;
    int ak;
    int al;
    int code;
};

struct cimproper_t {
    int code;
    real_t k;
    real_t phi0;
};

struct charge_t {
    int a;
    int code;
};

struct ccharge_t {
    int code;
    real_t charge;
};

struct atype_t {
    int a;
    int code;
};

struct catype_t {
    int code;
    real_t m;
    real_t aii_normal;
    real_t bii_normal;
    // real_t aii_polar;
    // real_t bii_polar;
    real_t aii_1_4;
    real_t bii_1_4;
};

struct vdw_pair_param_t {
    real_t a;
    real_t b;
};

struct topo_t {
    int solvent_type;
    real_t exclusion_radius;
    real_t solvent_radius;
    coord_t solute_center;
    coord_t solvent_center;
    real_t el14_scale;
    real_t coulomb_constant;
    int vdw_rule;  // 1=geometric, 2=arithmetic
};

struct cgrp_t {
    int n_atoms;
    int iswitch;
    int* a;
};

struct charge_group_t {
    int iswitch = 0;
    std::vector<int> atoms;
};

struct charge_group_config_t {
    int n_cgrps_solute = 0;
    int n_cgrps_solvent = 0;
    int iuse_switch_atom = 0;
    std::vector<charge_group_t> charge_groups;
};



struct q_angcouple_t { 
    int acode;
    int bcode;
}; // no use

struct q_cimproper_t {
    real_t k;
    real_t phi0;
}; // no use

struct q_elscale_t {
    int qi;
    int qj;
    real_t mu;
};

struct q_exclpair_t {
    int ai;
    int aj;
    int excl;
}; // no use

struct q_imprcouple_t {
    int icode;
    int bcode;
}; // no use

struct q_improper_t {
    int ai;
    int aj;
    int ak;
    int al;
    int code;
}; // no use

struct q_offdiag_t {
    int i;
    int j;
    int qk;
    int ql;
    real_t Aij;
    real_t muij;
}; // no use

struct q_shake_t {
    int ai;
    int aj;
    real_t dist;
}; // no use

struct q_softcore_t {
    real_t s;
}; // no use

struct q_softpair_t {
    int qi;
    int qj;
}; // no use

struct q_torcouple_t {
    int tcode;
    int bcode;
}; // no use

/* =============================================
 * == RESTRAINTS
 * =============================================
 */

struct restrseq_t {
    int ai;
    int aj;
    real_t k;
    bool ih;
    int to_center;  // Flag for restraining to geom. or mass center
};

struct restrpos_t {
    int a;
    int ipsi;
    coord_t x;
    coord_t k;
};

struct restrdis_t {
    int ai, aj;
    int ipsi;
    real_t d1, d2;
    real_t k;
    char itext[20], jtext[20];
};

struct restrang_t {
    int ai, aj, ak;
    int ipsi;
    real_t ang;
    real_t k;
};

struct restrwall_t {
    int ai, aj;
    real_t d, k, aMorse, dMorse;
    bool ih;
};

struct shell_t {
    int n_inshell;
    real_t theta_corr;
    real_t avtheta;
    real_t avn_inshell;
    real_t router;
    real_t dr;
    real_t cstb;
};

/* =============================================
 * == SHAKE
 * =============================================
 */

struct shake_bond_t {
    int ai;
    int aj;
    real_t dist2;
    bool ready;
};

/* =============================================
 * == CALCUTED IN THE INTEGRATION
 * =============================================
 */

struct vel_t {
    real_t x;
    real_t y;
    real_t z;
};

struct dvel_t {
    force_accum_t x;
    force_accum_t y;
    force_accum_t z;
};

struct E_bonded_t {
    real_t Ubond;
    real_t Uangle;
    real_t Utor;
    real_t Uimp;
};

struct E_nonbonded_t {
    real_t Ucoul;
    real_t Uvdw;
};

struct E_restraint_t {
    real_t Uradx;
    real_t Upolx;
    real_t Ufix;
    real_t Ushell;
    real_t Upres;
    real_t Urestr;
};

struct energy_t {
    real_t Ukin;
    real_t Upot;
    real_t Utot;
};
