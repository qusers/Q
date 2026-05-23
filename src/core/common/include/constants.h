// Boltzano's constant
#define Boltz 0.001986f

// Fortran max allowed line width, used in neighbor list
#define line_width 25

// Internally used time unit because of ??
#define time_unit 0.020462f

// Protein boundary force constant.
// TODO get force constant from md.inp
#define k_pshell 10.0

// Fixed proteins force constant.
#define k_fix 200.0

// Ratio of restrained protein shell that is free, rest is restrained. Has a default of 0.85
// TODO: get from md.inp
#define shell_default 0.85

// Definition of water shells
#define q_fortran_pi 3.1415927f
#define wpolr_layer 3.0001f
#define drouter 0.5f

// Number density of water in A measure
#define rho_water 0.0335f

// Once per how many steps theta_corr should be updated
#define itdis_update 100

// Shake convergence criterion (fraction of distance)
#define shake_tol 0.0001
#define shake_max_iter 1000

#define N_COLUMNS 10
#define COLUMN_WIDTH 25

#define VDW_GEOMETRIC 1
#define VDW_ARITHMETIC 2

// Thread block size
#define BLOCK_SIZE 8
