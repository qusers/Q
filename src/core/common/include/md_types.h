#pragma once

#include <string>
/* =============================================
 * == FROM MD FILE
 * =============================================
 */

struct md_t {
    // [MD]
    int steps;
    double stepsize;
    double temperature;
    char thermostat[40];
    double bath_coupling;
    int random_seed;
    double initial_temperature;
    bool shake_solvent;
    bool shake_solute;
    bool shake_hydrogens;
    bool lrf;
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
    double x;
    double y;
    double z;
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
    double charge;
};

struct atype_t {
    int a;
    int code;
};

struct catype_t {
    int code;
    double m;
    double aii_normal;
    double bii_normal;
    double aii_polar;
    double bii_polar;
    double aii_1_4;
    double bii_1_4;
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

struct q_angcouple_t {
    int acode;
    int bcode;
};

struct q_angle_t {
    int ai;
    int aj;
    int ak;
    int code;
};

struct q_atom_t {
    int a;
};

struct q_atype_t {
    int code;
};

struct q_bond_t {
    int ai;
    int aj;
    int code;
};

struct q_cangle_t {
    double kth;
    double th0;
};

struct q_catype_t {
    char name[10];
    double Ai;
    double Bi;
    double Ci;
    double ai;
    double Ai_14;
    double Bi_14;
    double m;
};

struct q_cbond_t {
    double kb;
    double b0;
};

struct q_charge_t {
    double q;
};

struct q_cimproper_t {
    double k;
    double phi0;
};

struct q_ctorsion_t {
    double k;
    double n;
    double d;
};

struct q_elscale_t {
    int qi;
    int qj;
    double mu;
};

struct q_exclpair_t {
    int ai;
    int aj;
    int excl;
};

struct q_imprcouple_t {
    int icode;
    int bcode;
};

struct q_improper_t {
    int ai;
    int aj;
    int ak;
    int al;
    int code;
};

struct q_offdiag_t {
    int i;
    int j;
    int qk;
    int ql;
    double Aij;
    double muij;
};

struct q_shake_t {
    int ai;
    int aj;
    double dist;
};

struct q_softcore_t {
    double s;
};

struct q_softpair_t {
    int qi;
    int qj;
};

struct q_torcouple_t {
    int tcode;
    int bcode;
};

struct q_torsion_t {
    int ai;
    int aj;
    int ak;
    int al;
    int code;
};

/* =============================================
 * == RESTRAINTS
 * =============================================
 */

struct restrseq_t {
    int ai;
    int aj;
    double k;
    bool ih;
    int to_center; // Flag for restraining to geom. or mass center
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

struct p_atom_t {
    int a;
};

struct vel_t {
    double x;
    double y;
    double z;
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
