#ifndef __SYSTEM_H__
#define __SYSTEM_H__

#define __PROFILING__
//#define DEBUG
#define VERBOSE

// Boltzano's constant
#define Boltz 0.001986

// Fortran max allowed line width, used in neighbor list
#define line_width 25

// Internally used time unit because of ??
#define time_unit 0.020462

// Protein boundary force constant.
// TODO get force constant from md.inp
#define k_pshell 10.0

// Fixed proteins force constant.
#define k_fix 200.0

// Ratio of restrained protein shell that is free, rest is restrained. Has a default of 0.85
// TODO: get from md.inp
#define shell_default 0.85

// Definition of water shells
#define wpolr_layer 3.0001
#define drouter 0.5

// Number density of water in A measure
#define rho_water 0.0335

// Once per how many steps theta_corr should be updated
#define itdis_update 100

// Shake convergence criterion (fraction of distance)
#define shake_tol 0.0001
#define shake_max_iter 1000

void init_variables();
void clean_variables();

/* =============================================
 * == DEVICE SETTINGS
 * =============================================
 */

// Thread block size
#define BLOCK_SIZE 8

/* =============================================
 * == GENERAL
 * =============================================
 */

extern int n_atoms;
extern int n_atoms_solute;
extern int n_patoms;
extern int n_qatoms;
extern int n_waters;
extern int n_molecules;

extern char base_folder[1024];

extern double dt;

extern bool run_gpu;

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
    double shell_radius; // Note: this is for the pshell
    double shell_force;  // Note: this is for the pshell
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

extern md_t md;
extern bool separate_scaling;

/* =============================================
 * == FROM TOPOLOGY FILE
 * =============================================
 */

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
};

struct cgrp_t {
    int n_atoms;
    int iswitch;
    int *a;
};

extern topo_t topo;

extern int n_angles;
extern int n_angles_solute;
extern int n_atypes;
extern int n_bonds;
extern int n_bonds_solute;
extern int n_cangles;
extern int n_catypes;
extern int n_cbonds;
extern int n_ccharges;
extern int n_charges;
extern int n_coords;
extern int n_cimpropers;
extern int n_ctorsions;
extern int n_excluded;
extern int n_impropers;
extern int n_impropers_solute;
extern int n_torsions;
extern int n_torsions_solute;
extern int n_excluded;
extern int n_cgrps_solute;
extern int n_cgrps_solvent;

extern angle_t *angles;
extern atype_t *atypes;
extern bond_t *bonds;
extern cangle_t *cangles;
extern catype_t *catypes;
extern cbond_t *cbonds;
extern charge_t *charges;
extern ccharge_t *ccharges;
extern cimproper_t *cimpropers;
extern ctorsion_t *ctorsions;
extern coord_t *coords_top;
extern coord_t *xcoords;
extern improper_t *impropers;
extern int *LJ_matrix;
extern torsion_t *torsions;
extern bool *excluded;
extern bool *heavy;
extern double *winv;
extern cgrp_t *charge_groups;

extern int *molecules;

/* =============================================
 * == FROM FEP FILE
 * =============================================
 */

extern int n_lambdas;
extern double *lambdas;

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

extern int n_qangcouples;
extern int n_qangles;
extern int n_qbonds;
extern int n_qcangles;
extern int n_qcatypes;
extern int n_qcbonds;
extern int n_qcimpropers;
extern int n_qctorsions;
extern int n_qelscales;
extern int n_qexclpairs;
extern int n_qimprcouples;
extern int n_qimpropers;
extern int n_qoffdiags;
extern int n_qshakes;
extern int n_qsoftcores;
extern int n_qsoftpairs;
extern int n_qtorcouples;
extern int n_qtorsions;

extern q_angcouple_t *q_angcouples;
extern q_atom_t *q_atoms;
extern q_cangle_t *q_cangles;
extern q_catype_t *q_catypes;
extern q_cbond_t *q_cbonds;
extern q_cimproper_t *q_cimpropers;
extern q_ctorsion_t *q_ctorsions;
extern q_offdiag_t *q_offdiags;
extern q_imprcouple_t *q_imprcouples;
extern q_softpair_t *q_softpairs;
extern q_torcouple_t *q_torcouples;

// NB. Arrays below are 2-dimensional!
extern q_angle_t *q_angles;
extern q_atype_t *q_atypes;
extern q_bond_t *q_bonds;
extern q_charge_t *q_charges;
extern q_elscale_t *q_elscales;
extern q_exclpair_t *q_exclpairs;
extern q_improper_t *q_impropers;
extern q_shake_t *q_shakes;
extern q_softcore_t *q_softcores;
extern q_torsion_t *q_torsions;

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

extern int n_restrseqs;
extern int n_restrspos;
extern int n_restrdists;
extern int n_restrangs;
extern int n_restrwalls;

extern restrseq_t *restrseqs;
extern restrpos_t *restrspos;
extern restrdis_t *restrdists;
extern restrang_t *restrangs;
extern restrwall_t *restrwalls;


// Protein shell layout. Defined once per run
extern bool *shell;

void init_pshells();
void init_restrseqs(char* filename);

struct shell_t {
    int n_inshell;
    double theta_corr;
    double avtheta;
    double avn_inshell;
    double router;
    double dr;
    double cstb;
};

// Total energy in the system. Defined once per run
extern double crgQtot;
extern double Dwmz, awmz;

// Water shell layout. Defined once per run
extern double *theta, *theta0, *tdum; //array size n_waters
extern int n_max_inshell, n_shells;
extern int **list_sh, **nsort; // array size (n_max_inshell, n_shells)
extern shell_t* wshells;

void init_water_sphere();
void init_wshells();

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

extern int n_shake_constraints, *mol_n_shakes;
extern shake_bond_t *shake_bonds;

void init_shake();

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

extern p_atom_t *p_atoms;
extern coord_t *coords;
extern vel_t* velocities;
extern dvel_t* dvelocities;
extern energy_t E_total;
extern energy_t *EQ_total;
extern E_bonded_t E_bond_p, E_bond_w, E_bond_q, *EQ_bond;
extern E_nonbonded_t E_nonbond_pp, E_nonbond_pw, E_nonbond_ww, E_nonbond_qx;
extern E_nonbonded_t *EQ_nonbond_qq, *EQ_nonbond_qp, *EQ_nonbond_qw, *EQ_nonbond_qx;
extern E_restraint_t E_restraint, *EQ_restraint;
extern double Temp, Tfree, Texcl;
extern double A_O, A_OO, B_O, B_OO, crg_ow, crg_hw; // TODO: don't keep this in system.cu?

void init_velocities();
void init_dvelocities();
void init_energies();

/* =============================================
 * == ENERGY & TEMPERATURE
 * =============================================
 */

extern double Ndegf, Ndegfree, Ndegf_solvent, Ndegf_solute, Ndegfree_solvent, Ndegfree_solute;
void calc_temperature();

/* =============================================
 * == INTEGRATION METHODS
 * =============================================
 */

void init_variables();
void clean_variables();
void write_header(const char *filename);
void write_coords(int iteration);
void write_velocities(int iteration);
void write_energies(int iteration);

void calc_integration();
void calc_integration_step(int iteration);

#endif /* __SYSTEM_H__ */