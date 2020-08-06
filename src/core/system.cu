#include "system.h"
#include "utils.h"
#include "parse.h"
#include "bonded.h"
#include "nonbonded.h"
#include "solvent.h"
#include "restraints.h"
#include "qatoms.h"

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <time.h>
#include <unistd.h>

/* =============================================
 * == GENERAL
 * =============================================
 */

int n_atoms;
int n_atoms_solute;
int n_patoms;
int n_qatoms;
int n_waters;

char base_folder[1024];
double dt, tau_T;

/* =============================================
 * == FROM MD FILE
 * =============================================
 */

md_t md;

/* =============================================
 * == FROM TOPOLOGY FILE
 * =============================================
 */

int n_coords;
int n_bonds;
int n_bonds_solute;
int n_cbonds;
int n_angles;
int n_angles_solute;
int n_cangles;
int n_torsions;
int n_torsions_solute;
int n_ctorsions;
int n_impropers;
int n_impropers_solute;
int n_cimpropers;
int n_charges;
int n_ccharges;
int n_atypes;
int n_catypes;
int n_ngbrs23;
int n_ngbrs14;

coord_t* coords_top;
bond_t* bonds;
cbond_t* cbonds;
angle_t* angles;
cangle_t* cangles;
torsion_t* torsions;
ctorsion_t* ctorsions;
improper_t* impropers;
cimproper_t* cimpropers;
charge_t* charges;
ccharge_t* ccharges;
atype_t* atypes;
catype_t* catypes;
ngbr23_t* ngbrs23;
ngbr14_t* ngbrs14;
bool *excluded;
bool *heavy;

topo_t topo;

/* =============================================
 * == FROM FEP FILE
 * =============================================
 */

int n_lambdas;
double *lambdas;

int n_qangcouples;
int n_qangles;
int n_qbonds;
int n_qcangles;
int n_qcatypes;
int n_qcbonds;
int n_qcimpropers;
int n_qctorsions;
int n_qelscales;
int n_qexclpairs;
int n_qimprcouples;
int n_qimpropers;
int n_qoffdiags;
int n_qshakes;
int n_qsoftpairs;
int n_qsoftcores;
int n_qtorcouples;
int n_qtorsions;

q_angcouple_t *q_angcouples;
q_atom_t *q_atoms;
q_cangle_t *q_cangles;
q_catype_t *q_catypes;
q_cbond_t *q_cbonds;
q_cimproper_t *q_cimpropers;
q_ctorsion_t *q_ctorsions;
q_offdiag_t *q_offdiags;
q_imprcouple_t *q_imprcouples;
q_softpair_t *q_softpairs;
q_torcouple_t *q_torcouples;

q_angle_t **q_angles;
q_atype_t **q_atypes;
q_bond_t **q_bonds;
q_charge_t **q_charges;
q_elscale_t **q_elscales;
q_exclpair_t **q_exclpairs;
q_improper_t **q_impropers;
q_shake_t **q_shakes;
q_softcore_t **q_softcores;
q_torsion_t **q_torsions;

// Remove bonds, angles, torsions and impropers which are excluded or changed in the FEP file
void shrink_topology() {
    int excluded;
    int ai = 0, bi = 0, ii = 0, ti = 0;
    int qai = 0, qbi = 0, qii = 0, qti = 0;

    excluded = 0;
    if (n_qangles > 0) {
        for (int i = 0; i < n_angles; i++) {
            if (angles[i].ai == q_angles[qai][0].ai
             && angles[i].aj == q_angles[qai][0].aj
             && angles[i].ak == q_angles[qai][0].ak) {
                qai++;
                excluded++;
            }
            else {
                angles[ai] = angles[i];
                ai++;
            }
        }
        n_angles -= excluded;
    }

    excluded = 0;
    if (n_qbonds > 0) {
        for (int i = 0; i < n_bonds; i++) {
            if (bonds[i].ai == q_bonds[qbi][0].ai
             && bonds[i].aj == q_bonds[qbi][0].aj) {
                qbi++;
                excluded++;
            }
            else {
                bonds[bi] = bonds[i];
                bi++;
            }
        }
        n_bonds -= excluded;
    }

    // excluded = 0;
    // for (int i = 0; i < n_impropers; i++) {
    //     if (impropers[i].ai == q_impropers[qai][0].ai
    //      && impropers[i].aj == q_impropers[qai][0].aj
    //      && impropers[i].ak == q_impropers[qai][0].ak
    //      && impropers[i].al == q_impropers[qai][0].al) {
    //         qii++;
    //         excluded++;
    //     }
    //     else {
    //         impropers[ii] = impropers[i];
    //         ii++;
    //     }
    // }
    // n_impropers -= excluded;

    excluded = 0;
    if (n_qtorsions > 0) {
        for (int i = 0; i < n_torsions; i++) {
            if (torsions[i].ai == q_torsions[qti][0].ai
             && torsions[i].aj == q_torsions[qti][0].aj
             && torsions[i].ak == q_torsions[qti][0].ak
             && torsions[i].al == q_torsions[qti][0].al) {
                qti++;
                excluded++;
            }
            else {
                torsions[ti] = torsions[i];
                ti++;
            }
        }
        n_torsions -= excluded;
    }

    // TODO: add exclusion pairs
}

/* =============================================
 * == CALCUTED IN THE INTEGRATION
 * =============================================
 */

coord_t* coords;
vel_t* velocities;
dvel_t* dvelocities;
energy_t energies;
double Temp = 0;
double Tscale = 1;

// Shake constrains
coord_t* xcoords;

// Water constants
double A_O = 0, A_OO = 0, B_O, B_OO, crg_ow, crg_hw, mu_w = 0;

void init_velocities() {
    velocities = (vel_t*) malloc(n_atoms * sizeof(vel_t));

    // If not previous value set, use a Maxwell distribution to fill velocities
    double kT = Boltz * md.initial_temperature;
    double sd, mass;
    for (int i = 0; i < n_atoms; i++) {
        mass = catypes[atypes[i].code - 1].m;
        sd = sqrt(kT / mass);

        velocities[i].x = gauss(0, sd);
        velocities[i].y = gauss(0, sd);
        velocities[i].z = gauss(0, sd);
    }
    
    //RANDOM SEED 1
    // velocities[0].x  = -1.9487889916926994E-002; velocities[0].y  = -3.9794110209903535E-004; velocities[0].z   = -7.5508872831782825E-004;
    // velocities[1].x  = 1.2857881903672132E-002; velocities[1].y  = -7.6644783763564926E-003; velocities[1].z   = -1.1772849886081924E-004;
    // velocities[2].x  = -1.6573999513626415E-002; velocities[2].y  = -4.9060190728011481E-003; velocities[2].z  = -3.1297180167266033E-003;
    // velocities[3].x  = 3.8335098604742421E-002;  velocities[3].y  = 1.7470716125507550E-002; velocities[3].z   = 1.7746920612499859E-002;
    // velocities[4].x  = -5.1831092590390600E-002;  velocities[4].y  = 3.3890986156632857E-002; velocities[4].z   = -1.8360181579421129E-002;
    // velocities[5].x  = -4.6344189702980876E-002; velocities[5].y  = 2.9483974897825213E-002; velocities[5].z   = 1.7345339668199650E-002;
    // velocities[6].x  = 3.3232237727582541E-002;  velocities[6].y  = 2.9621103153618156E-002; velocities[6].z   = 7.7704965685342450E-002;
    // velocities[7].x  = -6.5965875415921298E-003; velocities[7].y  = -1.3232621374638211E-002; velocities[7].z   = 4.6166681226082279E-002;
    // velocities[8].x  = -7.6695194596259841E-002;  velocities[8].y  = 1.3793159134795696E-002; velocities[8].z   = -7.5649865640426625E-003;
    // velocities[9].x  = 3.4339585945416480E-002;     velocities[9].y  = -1.9276301810268296E-002; velocities[9].z   = 2.1654891903647925E-003;
    // velocities[10].x = -8.6857064866742075E-003; velocities[10].y = 3.8524770102720978E-003; velocities[10].z  = 9.5000329534403842E-003;
    // velocities[11].x = 8.2276623425555515E-002;  velocities[11].y = -3.3500317868709793E-002; velocities[11].z  = -2.9493028790202912E-005;
    // velocities[12].x = 3.1561873458621555E-002; velocities[12].y = 1.0256001797406836E-002; velocities[12].z = 1.4081403694315784E-003;
    // velocities[13].x = 3.1497294075191563E-002; velocities[13].y = -3.2678383060448672E-002; velocities[13].z  = 6.1714353481433948E-002;
}

void init_dvelocities() {
    dvelocities = (dvel_t*) calloc(n_atoms, sizeof(dvel_t));
}

void init_xcoords() {
    xcoords = (coord_t*) malloc(n_atoms * sizeof(coord_t));
}

/* =============================================
 * == RESTRAINTS
 * =============================================
 */

// Array of length n_atoms
bool *shell;

int n_restrseqs;
int n_restrspos;
int n_restrdists;
int n_restrangs;
int n_restrwalls;

restrseq_t *restrseqs;
restrpos_t *restrspos;
restrdis_t *restrdists;
restrang_t *restrangs;
restrwall_t* restrwalls;

double crgQtot = 0;
double Dwmz, awmz;
 
 // Shell layout. Defined once per run
double *theta, *theta0, *tdum; //array size n_waters
int n_max_inshell, n_shells;
int **list_sh, **nsort; // array size (n_max_inshell, n_shells)
shell_t* wshells;


/* =============================================
 * == BOUNDARY RESTRAINTS
 * =============================================
 */

void init_water_sphere() {
    Dwmz = 0.26 * exp(-0.19 * (topo.solvent_radius - 15)) + 0.74;
    awmz = 0.2 / (1 + exp(0.4 * (topo.solvent_radius - 25))) + 0.3;

    printf("Dwmz = %f, awmz = %f\n", Dwmz, awmz);
}

//ONLY call if there are actually solvent atoms, or get segfaulted
void init_wshells() {
    int n_inshell;
    double drs, router, ri, dr, Vshell, rshell;
    if (mu_w == 0) {
        // Get water properties from first water molecule
        cbond_t cbondw = cbonds[bonds[n_atoms_solute].code-1];
        cangle_t canglew = cangles[angles[n_atoms_solute].code-1];

        ccharge_t ccharge_ow = ccharges[charges[n_atoms_solute].code - 1];
        crg_ow = ccharge_ow.charge;

        mu_w = -crg_ow * cbondw.b0 * cos(canglew.th0 / 2);
    }

    drs = wpolr_layer / drouter;

    n_shells = (int) floor(-0.5 + sqrt(2*drs + 0.25));
    wshells = (shell_t*) malloc(n_shells * sizeof(shell_t));

    printf("n_shells = %d\n", n_shells);

    router = topo.solvent_radius;
    n_max_inshell = 0;

    for (int i = 0; i < n_shells; i++) {
        wshells[i].avtheta = 0;
        wshells[i].avn_inshell = 0;
        wshells[i].router = router;
        dr = drouter * (i+1);
        ri = router - dr;
        wshells[i].dr = dr;
        Vshell = pow(router, 3) - pow(ri, 3);
        n_inshell = (int) floor(4 * M_PI / 3 * Vshell * rho_water);
        if (n_inshell > n_max_inshell) {
            n_max_inshell = n_inshell;
        }
        rshell = pow(0.5 * (pow(router, 3) + pow(ri, 3)), 1.0/3.0);

        // --- Note below: 0.98750 = (1-1/epsilon) for water
        wshells[i].cstb = crgQtot * 0.98750 / (rho_water * mu_w * 4 * M_PI * pow(rshell, 2));

        router -= dr;
    }

        // rc > wshells[n_shells-1].router - wshells[n_shells-1].dr
        printf("shell 0: (%f, %f). shell 1: (%f, %f). shell 2: (%f, %f).\n"
            , wshells[0].router, wshells[0].router - wshells[0].dr
            , wshells[1].router, wshells[1].router - wshells[1].dr
            , wshells[2].router, wshells[2].router - wshells[2].dr
        );

    n_max_inshell *= 1.5; // Make largest a little bigger just in case

    // Initialize arrays needed for bookkeeping
    theta = (double*) malloc(n_waters * sizeof(double));
    theta0 = (double*) malloc(n_waters * sizeof(double));
    tdum = (double*) malloc(n_waters * sizeof(double));

    list_sh = (int**) malloc(n_max_inshell * sizeof(int*));
    nsort = (int**) malloc(n_max_inshell * sizeof(int*));

    for (int i = 0; i < n_max_inshell; i++) {
        list_sh[i] = (int*) malloc(n_shells * sizeof(int));
        nsort[i] = (int*) malloc(n_shells * sizeof(int));
    }
}

void init_pshells() {
    double mass, r2, rin2;

    heavy = (bool*) malloc(n_atoms * sizeof(bool));
    shell = (bool*) malloc(n_atoms * sizeof(bool));
    rin2 = pow(shell_default * topo.exclusion_radius, 2);

    int n_heavy = 0, n_inshell = 0;

    for (int i = 0; i < n_atoms; i++) {
        mass = catypes[atypes[i].a-1].m;
        if (mass < 4.0) {
            heavy[i] = false;
        }
        else {
            heavy[i] = true;
            n_heavy++;
        }

        if (heavy[i] && !excluded[i] && i < n_atoms_solute) {
            r2 = pow(coords_top[i].x - topo.solute_center.x, 2) 
                + pow(coords_top[i].y - topo.solute_center.y, 2)
                + pow(coords_top[i].z - topo.solute_center.z, 2);
            if (r2 > rin2) {
                shell[i] = true;
                n_inshell++;
            }
            else {
                shell[i] = false;
            }
        }
    }

    printf("n_heavy = %d, n_inshell = %d\n", n_heavy, n_inshell);
}


void init_restrseqs() {
    n_restrseqs = 1;
    restrseqs = (restrseq_t*) malloc(1 * sizeof(restrseq_t));

    restrseq_t seq;
    seq.ai = 1;
    seq.aj = 14;
    seq.k = 1.0;
    seq.ih = 0;
    seq.to_center = 2;

    restrseqs[0] = seq;
}

/* =============================================
 * == CALCUTED IN THE INTEGRATION
 * =============================================
 */

p_atom_t *p_atoms;
energy_t *q_energies;

void init_patoms() {
    n_patoms = n_atoms_solute - n_qatoms;

    p_atoms = (p_atom_t*) malloc(n_patoms * sizeof(p_atom_t));

    // Loop through all solutes, adding a p atom to the list every time a non-q atom is encountered
    int pi = 0;
    int qi = 0;
    for (int i = 0; i < n_atoms_solute; i++) {
        if (n_qatoms > 0 && i == q_atoms[qi].a-1) {
            qi++;
        }
        else {
            p_atoms[pi].a = i+1;
            pi++;
        }
    }
}

/* =============================================
 * == ENERGY & TEMPERATURE
 * =============================================
 */

void calc_temperature() {
    double Ndegf = 3 * n_atoms;
    double Ekinmax = 1000.0 * Ndegf * Boltz * md.temperature / 2.0 / n_atoms;
    double ener;
    double mass_i;

    Temp = 0;
    for (int i = 0; i < n_atoms; i++) {
        mass_i = catypes[atypes[i].code - 1].m;
        ener = .5 * mass_i * (pow(velocities[i].x, 2) + pow(velocities[i].y, 2) + pow(velocities[i].z, 2));
        Temp += ener;
        if (ener > Ekinmax) {
            printf(">>> WARNING: hot atom %d: %f\n", i, ener/Boltz/3);
        }
    }

    energies.Ukin = Temp;

    Temp = 2.0 * Temp / Boltz / Ndegf;

    printf("Temp = %f\n", Temp);
    
    Tscale = sqrt(1 + (dt / tau_T) * (md.temperature / Temp - 1.0));
    // printf("Tscale = %f, tau_T = %f, Temp = %f\n", Tscale, tau_T, Temp);
}

/* =============================================
 * == INTEGRATION METHODS
 * =============================================
 */

void calc_leapfrog() {
    double mass_i, winv_i;
    for (int i = 0; i < n_atoms; i++) {
        mass_i = catypes[atypes[i].code - 1].m;

        winv_i = 1/mass_i;
        velocities[i].x = (velocities[i].x - dvelocities[i].x * dt * winv_i) * Tscale;
        velocities[i].y = (velocities[i].y - dvelocities[i].y * dt * winv_i) * Tscale;
        velocities[i].z = (velocities[i].z - dvelocities[i].z * dt * winv_i) * Tscale;

        // Prepare copy for shake
        xcoords[i].x = coords[i].x;
        xcoords[i].y = coords[i].y;
        xcoords[i].z = coords[i].z;

        coords[i].x += velocities[i].x * dt;
        coords[i].y += velocities[i].y * dt;
        coords[i].z += velocities[i].z * dt;
    }
}

// Write header (number of atoms) to output file
void write_header(const char *filename) {
    FILE * fp;

    char path[1024];
    sprintf(path, "%s/output/%s", base_folder, filename);

    fp = fopen(path, "w");
  
    fprintf(fp, "%d\n", n_atoms);
  
    fclose (fp);
}

// Write step number, coordinates of atoms to coordinate output file
void write_coords(int iteration) {
    if (iteration % md.trajectory != 0) return;
    FILE * fp;
    int i;

    char path[1024];
    sprintf(path, "%s/output/%s", base_folder, "coords.csv");

    fp = fopen(path, "a");

    fprintf(fp, "%d\n", iteration / md.trajectory);
    for(i = 0; i < n_atoms; i++) {
        fprintf(fp, "%f;%f;%f\n", coords[i].x, coords[i].y, coords[i].z);
    }
  
    fclose (fp);
}

// Write step number, velocities of atoms to coordinate output file
void write_velocities(int iteration) {
    if (iteration % md.trajectory != 0) return;
    FILE * fp;
    int i;

    char path[1024];
    sprintf(path, "%s/output/%s", base_folder, "velocities.csv");

    fp = fopen(path, "a");

    fprintf(fp, "%d\n", iteration / md.trajectory);
    for(i = 0; i < n_atoms; i++) {
        fprintf(fp, "%f;%f;%f\n", velocities[i].x, velocities[i].y, velocities[i].z);
    }
  
    fclose (fp);
}

// Write step number, energies of atoms to coordinate output file
void write_energies(int iteration) {
    if (iteration % md.energy != 0) return;
    FILE * fp;

    char path[1024];
    sprintf(path, "%s/output/%s", base_folder, "energies.csv");

    fp = fopen(path, "a");

    fprintf(fp, "%d\n", iteration / md.energy);

    for (int state = 0; state < n_lambdas; state++) {
        fprintf(fp, ">>> State %d\n", state);
        fprintf(fp, "Ubond = %f\n", q_energies[state].Ubond);
        fprintf(fp, "Uangle = %f\n", q_energies[state].Uangle);
        fprintf(fp, "Utor = %f\n", q_energies[state].Utor);
        fprintf(fp, "Ucoul = %f\n", q_energies[state].Ucoul);
        fprintf(fp, "Uvdw = %f\n", q_energies[state].Uvdw);
        fprintf(fp, "Urestr = %f\n", q_energies[state].Urestr);
        fprintf(fp, "Q-SUM = %f\n", q_energies[state].Upot);
    }

    fprintf(fp, ">>> Total\n");
    fprintf(fp, "Temp = %f\n", Temp);
    fprintf(fp, "Ubond = %f\n", energies.Ubond);
    fprintf(fp, "Uangle = %f\n", energies.Uangle);
    fprintf(fp, "Utor = %f\n", energies.Utor);
    fprintf(fp, "Uradx = %f\n", energies.Uradx);
    fprintf(fp, "Upolx = %f\n", energies.Upolx);
    fprintf(fp, "Ushell = %f\n", energies.Ushell);
    fprintf(fp, "Ufix = %f\n", energies.Ufix);
    fprintf(fp, "Upres = %f\n", energies.Upres);
    fprintf(fp, "Urestr = %f\n", energies.Urestr);
    fprintf(fp, "Ucoul = %f\n", energies.Ucoul);
    fprintf(fp, "Uvdw = %f\n", energies.Uvdw);
    fprintf(fp, "Ukin = %f\n", energies.Ukin);
    fprintf(fp, "Upot = %f\n", energies.Upot);
    fprintf(fp, "Utot = %f\n", energies.Utot);
  
    fclose (fp);
}

void calc_integration() {
    init_variables();

    for (int i = 0; i < md.steps; i++) {
        calc_integration_step(i);
    }
    
    clean_variables();
}

void calc_integration_step(int iteration) {
    printf("================================================\n");
    printf("== STEP %d\n", iteration);
    printf("================================================\n");

    // Reset derivatives
    for (int i = 0; i < n_atoms; i++) {
        dvelocities[i].x = 0;
        dvelocities[i].y = 0;
        dvelocities[i].z = 0;
    }
    energies.Ucoul = 0;
    energies.Uvdw = 0;
    for (int state = 0; state < n_lambdas; state++) {
        q_energies[state].Ucoul = 0;
        q_energies[state].Uvdw = 0;
    }

    // Determine temperature and kinetic energy
    calc_temperature();

    // Determine acceleration
    clock_t start = clock();

    // First solute interactions
    calc_angle_forces();
    calc_bond_forces();
    calc_torsion_forces();

    clock_t end_bonded = clock();

    calc_nonbonded_qp_forces();
    calc_nonbonded_pp_forces();

    clock_t end_nonbonded = clock();

    // Now solvent interactions
    if (n_waters > 0) {
        calc_nonbonded_ww_forces();
        calc_nonbonded_pw_forces();
        calc_nonbonded_qw_forces();
    }

    // Calculate restraints
    if (n_waters > 0) {
        calc_radix_w_forces();
        if (md.polarisation) {
            calc_polx_w_forces(iteration);
        }
    }
    calc_pshell_forces();
    calc_restrseq_forces();
    calc_restrdis_forces();

    // Q-Q nonbonded interactions
    calc_nonbonded_qq_forces();

    // Q-atom bonded interactions: loop over Q-atom states
    for (int state = 0; state < n_lambdas; state++) {
        calc_qangle_forces(state);
        calc_qbond_forces(state);
        calc_qtorsion_forces(state);
    }

    // Now apply leapfrog integration
    calc_leapfrog();

    // Recalculate temperature and kinetic energy for output
    calc_temperature();

    // Update energies with an average of all states
    for (int state = 0; state < n_lambdas; state++) {
        energies.Ubond += q_energies[state].Ubond * lambdas[state];
        energies.Uangle += q_energies[state].Uangle * lambdas[state];
        energies.Utor += q_energies[state].Utor * lambdas[state];
        energies.Ucoul += q_energies[state].Ucoul * lambdas[state];
        energies.Uvdw += q_energies[state].Uvdw * lambdas[state];

        // Update protein restraint energies with an average of all states
        energies.Upres += q_energies[state].Urestr * lambdas[state];

        // Update total Q-SUM for this state
        q_energies[state].Upot = q_energies[state].Ubond + q_energies[state].Uangle
            + q_energies[state].Utor + q_energies[state].Ucoul + q_energies[state].Uvdw
            + q_energies[state].Urestr;
    }

    // Update totals
    energies.Urestr = energies.Uradx + energies.Upolx + energies.Ushell + energies.Ufix + energies.Upres;
    energies.Upot = energies.Uangle + energies.Ubond + energies.Utor + energies.Ucoul + energies.Uvdw + energies.Urestr;
    energies.Utot = energies.Upot + energies.Ukin;  

    for (int state = 0; state < n_lambdas; state++) {
        printf(">>> State %d\n", state);
        printf("Ubond = %f\n", q_energies[state].Ubond);
        printf("Uangle = %f\n", q_energies[state].Uangle);
        printf("Utor = %f\n", q_energies[state].Utor);
        printf("Ucoul = %f\n", q_energies[state].Ucoul);
        printf("Uvdw = %f\n", q_energies[state].Uvdw);
        printf("Urestr = %f\n", q_energies[state].Urestr);
        printf("Q-SUM = %f\n", q_energies[state].Upot);
    }

    printf(">>> Total\n");
    printf("Tscale = %f\n", Tscale);
    printf("Ubond = %f\n", energies.Ubond);
    printf("Uangle = %f\n", energies.Uangle);
    printf("Utor = %f\n", energies.Utor);
    printf("Uradx = %f\n", energies.Uradx);
    printf("Upolx = %f\n", energies.Upolx);
    printf("Ushell = %f\n", energies.Ushell);
    printf("Ufix = %f\n", energies.Ufix);
    printf("Upres = %f\n", energies.Upres);
    printf("Urestr = %f\n", energies.Urestr);
    printf("Ucoul = %f\n", energies.Ucoul);
    printf("Uvdw = %f\n", energies.Uvdw);
    printf("Ukin = %f\n", energies.Ukin);
    printf("Upot = %f\n", energies.Upot);
    printf("Utot = %f\n", energies.Utot);


    // Profiler info
#ifdef __PROFILING__
    printf("Elapsed time for bonded forces: %f\n", (end_bonded-start) / (double)CLOCKS_PER_SEC );
    printf("Elapsed time for non-bonded forces: %f\n", (end_nonbonded-end_bonded) / (double)CLOCKS_PER_SEC);
#endif /* __PROFILING__ */

    // Append output files
    write_coords(iteration);
    write_velocities(iteration);
    write_energies(iteration);
}

void init_variables() {
    // From MD file
    init_md("md.csv");

    dt = time_unit * md.stepsize;
    tau_T = time_unit * md.bath_coupling;

    // From topology file
    init_topo("topo.csv");
    
    init_angles("angles.csv");
    init_atypes("atypes.csv");
    init_bonds("bonds.csv");
    init_cangles("cangles.csv");
    init_catypes("catypes.csv");
    init_cbonds("cbonds.csv");
    init_ccharges("ccharges.csv");
    init_charges("charges.csv");
    init_cimpropers("cimpropers.csv");
    init_coords("coords.csv");
    init_ctorsions("ctorsions.csv");
    init_excluded("excluded.csv");
    init_impropers("impropers.csv");
    init_torsions("torsions.csv");
    init_ngbrs14("ngbrs14.csv");
    init_ngbrs23("ngbrs23.csv");
    // init_restrseqs();

    // From FEP file
    init_qangcouples("q_angcouples.csv");
    init_qatoms("q_atoms.csv");
    init_qcangles("q_cangles.csv");
    init_qcatypes("q_catypes.csv");
    init_qcbonds("q_cbonds.csv");
    init_qcimpropers("q_cimpropers.csv");
    init_qctorsions("q_ctorsions.csv");
    init_qoffdiags("q_offdiags.csv");
    init_qimprcouples("q_imprcouples.csv");
    init_qsoftpairs("q_softpairs.csv");
    init_qtorcouples("q_torcouples.csv");

    init_qangles("q_angles.csv");
    init_qatypes("q_atypes.csv");
    init_qbonds("q_bonds.csv");
    init_qcharges("q_charges.csv");
    init_qelscales("q_elscales.csv");
    init_qexclpairs("q_exclpairs.csv");
    init_qimpropers("q_impropers.csv");
    init_qshakes("q_shakes.csv");
    init_qsoftcores("q_softcores.csv");
    init_qtorsions("q_torsions.csv");

    shrink_topology();

    // Init random seed from MD file
    srand(md.random_seed);

    // From calculation in the integration
    init_patoms();
    init_velocities();
    init_dvelocities();
    init_xcoords();

    // From input file
    // init_icoords("i_coords.csv");
    // init_ivelocities("i_velocities.csv");    

    // Init waters, boundary restrains
    n_waters = (n_atoms - n_atoms_solute) / 3;
    if (n_waters > 0) {
        init_water_sphere();
        init_wshells();
    }
    init_pshells();

    // Init energy
    energies.Ubond = 0;
    energies.Uangle = 0;
    energies.Utor = 0;
    energies.Ucoul = 0;
    energies.Uvdw = 0;
    energies.Ukin = 0;
    energies.Upot = 0;
    energies.Uradx = 0;
    energies.Upolx = 0;
    energies.Ushell = 0;
    energies.Ufix = 0;
    energies.Urestr = 0;

    q_energies = (energy_t*) calloc(n_lambdas, sizeof(energy_t));

    // Write header to file
    write_header("coords.csv");
    write_header("velocities.csv");
    write_header("energies.csv");
}

void clean_variables() {
    // From topology file
    free(angles);
    free(atypes);
    free(bonds);
    free(cangles);
    free(catypes);
    free(cbonds);
    free(ccharges);
    free(charges);
    free(cimpropers);
    free(coords);
    free(coords_top);
    free(ctorsions);
    free(excluded);
    free(heavy);
    free(impropers);
    free(torsions);
    free(ngbrs14);
    free(ngbrs23);

    // From FEP file
    free(q_angcouples);
    free(q_atoms);
    free(q_cangles);
    free(q_catypes);
    free(q_cbonds);
    free(q_cimpropers);
    free(q_ctorsions);
    free(q_offdiags);
    free(q_imprcouples);
    free(q_softpairs);
    free(q_torcouples);
    for (int i = 0; i < n_qatoms; i++) {
        free(q_atypes[i]);
        free(q_charges[i]);
    }
    free(q_atypes);
    free(q_charges);
    for (int i = 0; i < n_qangles; i++) {
        free(q_angles[i]);
    }
    free(q_angles);
    for (int i = 0; i < n_qbonds; i++) {
        free(q_bonds[i]);
    }
    free(q_bonds);
    for (int i = 0; i < n_qelscales; i++) {
        free(q_elscales[i]);
    }
    free(q_elscales);
    for (int i = 0; i < n_qexclpairs; i++) {
        free(q_exclpairs[i]);
    }
    free(q_exclpairs);
    for (int i = 0; i < n_qimpropers; i++) {
        free(q_impropers[i]);
    }
    free(q_impropers);
    for (int i = 0; i < n_qshakes; i++) {
        free(q_shakes[i]);
    }
    free(q_shakes);
    for (int i = 0; i < n_qsoftcores; i++) {
        free(q_softcores[i]);
    }
    free(q_softcores);
    for (int i = 0; i < n_qtorsions; i++) {
        free(q_torsions[i]);
    }
    free(q_torsions);

    // Restraints
    if (n_waters > 0) {
        free(wshells);
        free(theta);
        free(theta0);
        free(tdum);
    
        for (int i = 0; i < n_max_inshell; i++) {
            free(list_sh[i]);
            free(nsort[i]);
        }

        free(list_sh);
        free(nsort);
    }
    free(restrseqs);
    free(shell);

    // From calculation in the integration
    free(velocities);
    free(dvelocities);
    free(xcoords);

    // Energies & temperature
    free(q_energies);
}