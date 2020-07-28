#ifndef __SYSTEM_H__
#define __SYSTEM_H__

//#define __PROFILING__

// Coulomb's constant, TODO get this from topology file
#define Coul 332.0716

// Boltzano's constant
#define Boltz 0.001986

// Maxwell temperature, TODO get this from md.inp
#define Tmaxw 300.0

// Target temperature, TODO get this from md.inp
#define Temp0 300.0

// Fortran max allowed line width, used in neighbor list
#define line_width 25

// Timestep definitions, 1 step = 10E-15s
// TODO get step_size, steps from md.inp
//#define step_size_unit 10E-5
#define step_size 1
//#define steps 10

// Internally used time unit because of ??
#define time_unit 0.020462

// Coupling constant used in temperature calculation
// TODO get bath_coupling from md.inp
#define c_bath_coupling 1.0

// Surface inward harmonic force constant. Has a default of 60
// TODO get force constant from md.inp
#define k_wsphere 60.0

// Polarization force constant. Has a default of 20
// TODO get force constant from md.inp
#define k_wpol 20.0

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

void init_variables();
void clean_variables();

/* =============================================
 * == GENERAL
 * =============================================
 */

extern int n_atoms;
extern int n_atoms_solute;
extern int n_waters;

extern char base_folder[1024];

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
    bool shake_hydrogens;
    bool lrf;
    // [cut-offs]
    double solute_solute;
    double solvent_solvent;
    double solute_solvent;
    double q_atom;
    // [sphere]
    double shell_radius;
    double shell_force;
    // [solvent]
    double radial_force;
    bool polarization;
    double polarization_force;
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
    double n;
    double d;
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

struct ngbr23_t {
    int ai;
    int aj;
};

struct ngbr14_t {
    int ai;
    int aj;
};

struct topo_t {
    int solvent_type;
    double exclusion_radius;
    double solvent_radius;
    coord_t solute_center;
    coord_t solvent_center;
};

extern topo_t topo;

extern int n_angles;
extern int n_angles_solute;
extern int n_atypes;
extern int n_bonds;
extern int n_bonds_solute;
extern int n_cangles;
extern int n_cbonds;
extern int n_ccharges;
extern int n_charges;
extern int n_coords;
extern int n_cimpropers;
extern int n_ctorsions;
extern int n_impropers;
extern int n_impropers_solute;
extern int n_ngbrs14;
extern int n_ngbrs23;
extern int n_torsions;
extern int n_torsions_solute;

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
extern improper_t *impropers;
extern ngbr14_t *ngbrs14;
extern ngbr23_t *ngbrs23;
extern torsion_t *torsions;
extern bool *excluded;
extern bool *heavy;

/* =============================================
 * == Q ATOMS
 * =============================================
 */

extern int n_lambdas;
extern double *lambdas;

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
 * == CALCUTED IN THE INTEGRATION
 * =============================================
 */

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

struct energy_t {
    double Ubond;
    double Uangle;
    double Utor;
    double Ucoul;
    double Uvdw;
    double Ukin;
    double Upot;
    double Utot;
    double Uradx;
    double Upolx;
    double Ufix;
    double Ushell;
    double Upres;
    double Urestr;
};

extern coord_t *coords;
extern dvel_t* dvelocities;
extern energy_t energies;
extern double Temp;
extern double A_OO, B_OO, crg_ow, crg_hw; // TODO: don't keep this in system.cu?

void init_velocities();
void init_dvelocities();
void init_energies();

/* =============================================
 * == INTEGRATION METHODS
 * =============================================
 */

void calc_integration();
void calc_integration_step(int iteration);

#endif /* __SYSTEM_H__ */