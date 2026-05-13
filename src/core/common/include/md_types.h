#pragma once

#include <string>
#include <vector>

#include "common/include/precision.h"

struct md_t {
    // [MD]
    int steps;
    double stepsize;
    double temperature;
    std::string thermostat;
    double bath_coupling;
    int random_seed;
    double initial_temperature;
    bool shake_solvent;
    bool shake_solute;
    bool shake_hydrogens;
    bool lrf;
    bool separate_scaling = false;
    bool charge_groups;
    // [cut-offs]
    double solute_solute;
    double solvent_solvent;
    double solute_solvent;
    double q_atom;
    // [sphere]
    double shell_radius;  // Note: this is for the pshell
    double shell_force;   // Note: this is for the pshell
    // [solvent]
    double radial_force;
    bool polarisation;
    double polarisation_force;
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
    double kb;
    double b0;
};

struct angle_t {
    int ai;
    int aj;
    int ak;
    int code;
};

struct cangle_t {
    int code;
    double kth;
    double th0;
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
    double k;
    double n;
    double d;
    double paths;
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
    double k;
    double phi0;
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
    double m;
    real_t aii_normal;
    real_t bii_normal;
    // double aii_polar;
    // double bii_polar;
    real_t aii_1_4;
    real_t bii_1_4;
};

struct vdw_pair_param_t {
    real_t a;
    real_t b;
};

struct topo_t {
    int solvent_type;
    double exclusion_radius;
    double solvent_radius;
    coord_t solute_center;
    coord_t solvent_center;
    double el14_scale;
    double coulomb_constant;
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
    double k;
    double phi0;
}; // no use

struct q_elscale_t {
    int qi;
    int qj;
    double mu;
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
    double Aij;
    double muij;
}; // no use

struct q_shake_t {
    int ai;
    int aj;
    double dist;
}; // no use

struct q_softcore_t {
    double s;
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
    double k;
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
    double d1, d2;
    double k;
    char itext[20], jtext[20];
};

struct restrang_t {
    int ai, aj, ak;
    int ipsi;
    double ang;
    double k;
};

struct restrwall_t {
    int ai, aj;
    double d, k, aMorse, dMorse;
    bool ih;
};

struct shell_t {
    int n_inshell;
    double theta_corr;
    double avtheta;
    double avn_inshell;
    double router;
    double dr;
    double cstb;
};

/* =============================================
 * == SHAKE
 * =============================================
 */

struct shake_bond_t {
    int ai;
    int aj;
    double dist2;
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
    double x;
    double y;
    double z;
};

struct E_bonded_t {
    double Ubond;
    double Uangle;
    double Utor;
    double Uimp;
};

struct E_nonbonded_t {
    double Ucoul;
    double Uvdw;
};

struct E_restraint_t {
    double Uradx;
    double Upolx;
    double Ufix;
    double Ushell;
    double Upres;
    double Urestr;
};

struct energy_t {
    double Ukin;
    double Upot;
    double Utot;
};
