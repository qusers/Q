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
    // todo: will remove later
    bool shake_solvent;
    bool shake_solute;
    bool shake_hydrogens;

    bool constraint_solvent;
    bool constraint_solute;
    bool constraint_hydrogens;
    std::string solute_constraint_algorithm;
    std::string solvent_constraint_algorithm;

    bool lrf;
    bool separate_scaling = false;
    bool charge_groups;

    bool hmr = false;
    double hmr_target_mass = 3.024;

    // [cut-offs]
    double solute_solute;
    double solvent_solvent;
    double solute_solvent;
    double q_atom;
    // [sphere]
    bool exclude_bonded = false;
    double shell_radius;  // Note: this is for the pshell
    double shell_force;   // Note: this is for the pshell
    // [solvent]
    double radial_force;
    bool polarisation;
    bool charge_correction;
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

template <class T>
struct Real3 {
    T x, y, z;
};

template <class T>
struct Real2 {
    T x, y;
};

using coord_t = Real3<double>;
using real_t3 = Real3<real_t>;
using real_t2 = Real2<real_t>;

struct bond_t {
    int ai;
    int aj;
    int code;
};

struct cbond_t {
    int code = 0;
    double kb = 0.0;
    double b0 = 0.0;

    double D = 0.0;
    double a = 0.0;
    bool is_morse = false;
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
    // double aii_polar;
    // double bii_polar;
    double aii_1_4;
    double bii_1_4;

    double soft_C = 0.0;
    double soft_a = 0.0;
};

struct vdw_atom_param_t {
    real_t aii_normal;
    real_t bii_normal;
    real_t aii_14;
    real_t bii_14;
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
    int angle_index = 0;  // 1-based row index in [change_angles]
    int bond_index = 0;   // 1-based row index in [change_bonds]
};  // no use

struct q_cimproper_t {
    double k;
    double phi0;
};  // no use

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
    int improper_index = 0;  // 1-based row in [change_impropers]
    int bond_index = 0;      // 1-based row in [change_bonds]
    bool inverted = false;   // true: use 1 - gamma
};

struct q_improper_t {
    int ai;
    int aj;
    int ak;
    int al;
    int code;
};  // no use

struct q_offdiag_t {
    int i = 0;   // 1-based state index
    int j = 0;   // 1-based state index
    int qk = 0;  // 1-based Q-atom index in [atoms]
    int ql = 0;  // 1-based Q-atom index in [atoms]
    double Aij = 0.0;
    double muij = 0.0;
};

struct offdiag_term_t {
    int state_i = 0;  // 1-based, retained for .en output
    int state_j = 0;
    int ak = -1;  // 0-based physical atom index
    int al = -1;
    double A = 0.0;
    double mu = 0.0;
};

struct offdiag_record_t {
    int state_i = 0;  // 1-based
    int state_j = 0;
    double Hij = 0.0;
    double rkl = 0.0;
};

struct q_shake_t {
    int ai;
    int aj;
    double dist;
};  // no use

struct q_softcore_t {
    double s;
};  // no use

struct q_softpair_t {
    int qi;  // 1-based Q-atom index
    int qj;  // 1-based Q-atom index
};

struct q_torcouple_t {
    int torsion_index = 0;  // 1-based row index in [change_torsions]
    int bond_index = 0;     // 1-based row index in [change_bonds]
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

/* =============================================
 * == CALCUTED IN THE INTEGRATION
 * =============================================
 */

struct vel_t {
    double x;
    double y;
    double z;
};

struct dvel_t {
    force_accum_t x;
    force_accum_t y;
    force_accum_t z;
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

enum class SpecialPairKind : int {
    SoftRepulsion,
    LJ14,  // Reserved; not enabled yet
};

struct special_pair_term_t {
    int ai = -1;  // 0-bas
    int aj = -1;
    int state = -1;  // -1: weight 1; otherwise lambdas[state]

    SpecialPairKind kind = SpecialPairKind::SoftRepulsion;

    int coul_slot = -1;
    int vdw_slot = -1;

    double coulomb_coeff = 0.0;

    // SoftRepulsion: p1 = Ci * Cj, p2 = ai * aj
    // LJ14: p1 = combined A14, p2 = combined B14
    double p1 = 0.0;
    double p2 = 0.0;
};
