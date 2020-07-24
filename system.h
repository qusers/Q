#ifndef __SYSTEM_H__
#define __SYSTEM_H__

// Coulomb's constant, TODO get this from topology file
#define Coul 332.0716

// Boltzano's constant
#define Boltz 0.001986

// Maxwell temperature, TODO get this from md.inp
#define Tmaxw 300.0

// Initial temperature, TODO get this from md.inp
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
#define bath_coupling 1.0

// Center of the water sphere
// TODO get sphere coordinates from topology
#define centerX 0.0
#define centerY 0.0
#define centerZ 0.0
#define rwater 15.000

// Surface inward harmonic force constant. Has a default of 60
// TODO get force constant from md.inp
#define k_wsphere 60.0

// Polarization force constant. Has a default of 20
// TODO get force constant from md.inp
#define k_wpol 20.0

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

extern int n_angles;
extern int n_bonds;
extern int n_ngbrs14;
extern int n_ngbrs23;
extern int n_torsions;

extern angle_t* angles;
extern atype_t* atypes;
extern bond_t* bonds;
extern cangle_t* cangles;
extern catype_t* catypes;
extern cbond_t* cbonds;
extern charge_t* charges;
extern ccharge_t* ccharges;
extern ctorsion_t* ctorsions;
extern coord_t* coords;
extern ngbr14_t* ngbrs14;
extern ngbr23_t* ngbrs23;
extern torsion_t* torsions;

void init_coords(char* filename);
void init_bonds(char* filename);
void init_cbonds(char* filename);
void init_angles(char* filename);
void init_cangles(char* filename);
void init_torsions(char* filename);
void init_ctorsions(char* filename);
void init_impropers(char* filename);
void init_cimpropers(char* filename);
void init_charges(char* filename);
void init_ccharges(char* filename);
void init_ngbr14s(char* filename);
void init_ngbr23s(char* filename);
void init_catypes(char* filename);
void init_atypes(char* filename);

/* =============================================
 * == RESTRAINTS
 * =============================================
 */

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

// Shell layout. Defined once per run
extern double *theta, *theta0, *tdum; //array size n_waters
extern int n_max_inshell, n_shells;
extern int **list_sh, **nsort; // array size (n_max_inshell, n_shells)
extern shell_t* wshells;

void init_water_sphere();
void init_shells();

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
    double Urestr;
};

extern dvel_t* dvelocities;
extern energy_t energies;
extern double Temp;
extern double A_OO, B_OO, crg_ow, crg_hw; // TODO: don't keep this in system.cu?

void init_velocities();
void init_dvelocities();
void init_energies();
void init_water_sphere();

/* =============================================
 * == INTEGRATION METHODS
 * =============================================
 */

void calc_integration(int iteration);

#endif /* __SYSTEM_H__ */