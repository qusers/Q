#include "system.h"
#include "utils.h"
#include "parse.h"
#include "bonded.h"
#include "nonbonded.h"
#include "solvent.h"
#include "restraints.h"
#include "qatoms.h"
#include "shake.h"

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
int n_molecules = 0;

char base_folder[1024];
double dt, tau_T;

bool run_gpu = false;

/* =============================================
 * == FROM MD FILE
 * =============================================
 */

md_t md;
bool separate_scaling = false;

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
int n_excluded;
int n_cgrps_solute;
int n_cgrps_solvent;

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
int *LJ_matrix;
bool *excluded;
bool *heavy;
int *molecules;
double *winv;
cgrp_t *charge_groups;

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

q_angle_t *q_angles;
q_atype_t *q_atypes;
q_bond_t *q_bonds;
q_charge_t *q_charges;
q_elscale_t *q_elscales;
q_exclpair_t *q_exclpairs;
q_improper_t *q_impropers;
q_shake_t *q_shakes;
q_softcore_t *q_softcores;
q_torsion_t *q_torsions;

// Remove bonds, angles, torsions and impropers which are excluded or changed in the FEP file
void exclude_qatom_definitions() {
    int excluded;
    int ai = 0, bi = 0, ii = 0, ti = 0;
    int qai = 0, qbi = 0, qii = 0, qti = 0;

    excluded = 0;
    if (n_qangles > 0) {
        for (int i = 0; i < n_angles; i++) {
            if (angles[i].ai == q_angles[qai].ai
             && angles[i].aj == q_angles[qai].aj
             && angles[i].ak == q_angles[qai].ak) {
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
            if (bonds[i].ai == q_bonds[qbi].ai
             && bonds[i].aj == q_bonds[qbi].aj) {
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
            if (torsions[i].ai == q_torsions[qti].ai
             && torsions[i].aj == q_torsions[qti].aj
             && torsions[i].ak == q_torsions[qti].ak
             && torsions[i].al == q_torsions[qti].al) {
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

void exclude_all_atoms_excluded_definitions() {
    int n_excl;
    int ai = 0, bi = 0, ii = 0, ti = 0;

    // n_excl = 0;
    // for (int i = 0; i < n_angles; i++) {
    //     if (excluded[angles[i].ai - 1]
    //         && excluded[angles[i].aj - 1]
    //         && excluded[angles[i].ak - 1]) {
    //         n_excl++;
    //     }
    //     else {
    //         angles[ai] = angles[i];
    //         ai++;
    //     }
    // }
    // printf("original: %d. # excluded angles: %d\n", n_angles, n_excl);
    // n_angles -= n_excl;

    // n_excl = 0;
    // for (int i = 0; i < n_bonds; i++) {
    //     if (excluded[bonds[i].ai - 1]
    //         && excluded[bonds[i].aj - 1]) {
    //         n_excl++;
    //     }
    //     else {
    //         bonds[bi] = bonds[i];
    //         bi++;
    //     }
    // }
    // printf("original: %d. # excluded bonds: %d\n", n_bonds, n_excl);
    // n_bonds -= n_excl;

    n_excl = 0;
    for (int i = 0; i < n_impropers; i++) {
        if (excluded[impropers[i].ai - 1]
            && excluded[impropers[i].aj - 1]
            && excluded[impropers[i].ak - 1]
            && excluded[impropers[i].al - 1]) {
            n_excl++;
        }
        else {
            impropers[ii] = impropers[i];
            ii++;
        }
    }
    printf("original: %d. # excluded impropers: %d\n", n_impropers, n_excl);
    n_impropers -= n_excl;

    n_excl = 0;
    for (int i = 0; i < n_torsions; i++) {
        if (excluded[torsions[i].ai - 1]
            && excluded[torsions[i].aj - 1]
            && excluded[torsions[i].ak - 1]
            && excluded[torsions[i].al - 1]) {
            n_excl++;
        }
        else {
            torsions[ti] = torsions[i];
            ti++;
        }
    }
    printf("original: %d. # excluded torsions: %d\n", n_torsions, n_excl);
    n_torsions -= n_excl;
}

void exclude_shaken_definitions() {
    int excluded;
    int bi = 0;
    int si = 0;
    int ang_i = 0;
    int ai, aj;

    excluded = 0;
    if (n_shake_constraints > 0) {
        for (int i = 0; i < n_bonds; i++) {
            if (bonds[i].ai == shake_bonds[si].ai
             && bonds[i].aj == shake_bonds[si].aj) {
                si++;
                excluded++;
            }
            else {
                bonds[bi] = bonds[i];
                bi++;
            }
        }
        n_bonds -= excluded;
    }

    excluded = 0;
    if (n_shake_constraints > 0) {
        for (int i = 0; i < n_shake_constraints; i++) {
            ai = shake_bonds[i].ai;
            aj = shake_bonds[i].aj;
            for (int j = 0; j < n_angles; j++) {
                if ( (angles[j].ai == ai && angles[j].aj == aj)
                    || (angles[j].ai == aj && angles[j].ak == aj) ) {
                    angles[j].code = 0;
                    break;
                }
            }
        }

        for (int i = 0; i < n_angles; i++) {
            if (angles[i].code == 0) {
                excluded++;
            }
            else {
                angles[ang_i] = angles[i];
                ang_i++;
            }
        }
    }
}

/* =============================================
 * == CALCUTED IN THE INTEGRATION
 * =============================================
 */

coord_t* coords;
vel_t* velocities;
dvel_t* dvelocities;
double Temp = 0;
double Texcl = 0;
double Tfree = 0;
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
}

void init_dvelocities() {
    dvelocities = (dvel_t*) calloc(n_atoms, sizeof(dvel_t));
}

void init_xcoords() {
    xcoords = (coord_t*) malloc(n_atoms * sizeof(coord_t));
}

void init_inv_mass() {
    winv = (double*) malloc(n_atoms * sizeof(double));
    for (int ai = 0; ai < n_atoms; ai++) {
        winv[ai] = 1 / catypes[atypes[ai].code-1].m;
    }
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

    n_max_inshell = n_waters; // Make largest a little bigger just in case

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

    heavy = (bool*) calloc(n_atoms, sizeof(bool));
    shell = (bool*) calloc(n_atoms, sizeof(bool));
    rin2 = pow(shell_default * topo.exclusion_radius, 2);

    int n_heavy = 0, n_inshell = 0;

    for (int i = 0; i < n_atoms; i++) {
        mass = catypes[atypes[i].code-1].m;
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

void init_pshells_with_charge_groups() {
    double mass, r2, rin2;

    heavy = (bool*) calloc(n_atoms, sizeof(bool));
    shell = (bool*) calloc(n_atoms, sizeof(bool));
    rin2 = pow(md.shell_radius, 2);

    int n_heavy = 0, n_inshell = 0;

    for (int i = 0; i < n_atoms; i++) {
        mass = catypes[atypes[i].code-1].m;
        if (mass < 4.0) {
            heavy[i] = false;
        }
        else {
            heavy[i] = true;
            n_heavy++;
        }
    }

    for (int grp = 0; grp < n_cgrps_solute; grp++) {
        cgrp_t cgrp = charge_groups[grp];
        int i = cgrp.iswitch-1;
        if (heavy[i] && !excluded[i] && i < n_atoms_solute) {
            r2 = pow(coords_top[i].x - topo.solute_center.x, 2) 
                + pow(coords_top[i].y - topo.solute_center.y, 2)
                + pow(coords_top[i].z - topo.solute_center.z, 2);
            bool switch_atom_in_shell = r2 > rin2;
            for (int j = 0; j < cgrp.n_atoms; j++) {
                shell[cgrp.a[j]-1] = switch_atom_in_shell;
                if (switch_atom_in_shell) {
                    n_inshell++;
                }
            }
        }
    }

    printf("(with charge groups): n_heavy = %d, n_inshell = %d\n", n_heavy, n_inshell);
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
 * == SHAKE
 * =============================================
 */

int n_shake_constraints, *mol_n_shakes;
shake_bond_t *shake_bonds;
 
void init_shake() {
    int ai, aj;
    int mol = 0;
    int shake;
    int n_solute_shake_constraints = 0;
    double excl_shake = 0;

    n_shake_constraints = 0;
    mol_n_shakes = (int*) calloc(n_molecules, sizeof(int));
    
    for (int bi = 0; bi < n_bonds; bi++) {
        ai = bonds[bi].ai-1;
        aj = bonds[bi].aj-1;

        while(mol+1 < n_molecules && ai+1 >= molecules[mol+1]) {
            // new molecule
            mol += 1;
        }

        if ( (md.shake_hydrogens && (!heavy[ai] || !heavy[aj]))
            || (md.shake_solute && ai+1 <= n_atoms_solute) 
            || (md.shake_solvent && ai+1 > n_atoms_solute) ) {
            mol_n_shakes[mol]++;
            n_shake_constraints++;

            if (excluded[ai]) excl_shake += 0.5;
            if (excluded[aj]) excl_shake += 0.5;
    
        }
    }

    shake_bonds = (shake_bond_t*) malloc(n_shake_constraints * sizeof(shake_bond_t));
    mol = 0;
    shake = 0;
    for (int bi = 0; bi < n_bonds; bi++) {
        ai = bonds[bi].ai-1;
        aj = bonds[bi].aj-1;

        while(ai+1 >= molecules[mol+1]) {
            // new molecule
            mol += 1;
        }

        if ( (md.shake_hydrogens && (!heavy[ai] || !heavy[aj]))
            || (md.shake_solute && ai+1 <= n_atoms_solute) 
            || (md.shake_solvent && ai+1 > n_atoms_solute) ) {
            shake_bonds[shake].ai = ai+1;
            shake_bonds[shake].aj = aj+1;
            shake_bonds[shake].dist2 = pow(cbonds[bonds[bi].code-1].b0, 2);
            shake++;
        }
    }

    // Get total number of shake constraints in solute (used for separate scaling of temperatures)
    for (int i = 0; i < n_molecules - n_waters; i++) {
        n_solute_shake_constraints += mol_n_shakes[i];
    }

    Ndegf = 3 * n_atoms - n_shake_constraints;
    Ndegfree = Ndegf - 3 * n_excluded + excl_shake;

    Ndegf_solvent = Ndegf - 3 * n_atoms_solute  + n_solute_shake_constraints;
    Ndegf_solute = Ndegf - Ndegf_solvent;

    Ndegfree_solvent = Ndegfree - (n_shake_constraints - n_solute_shake_constraints);
    Ndegfree_solute = Ndegfree - Ndegfree_solvent;

    printf("n_shake_constrains = %d, n_solute_shake_constraints = %d, excl_shake = %f\n", n_shake_constraints, n_solute_shake_constraints, excl_shake);

    if (Ndegfree_solvent * Ndegfree_solute == 0) {
        separate_scaling = false;
    }
    else {
        separate_scaling = true;
    }
}

/* =============================================
 * == CALCUTED IN THE INTEGRATION
 * =============================================
 */

p_atom_t *p_atoms;
energy_t E_total;
energy_t *EQ_total;
E_bonded_t E_bond_p, E_bond_w, E_bond_q, *EQ_bond;
E_nonbonded_t E_nonbond_pp, E_nonbond_pw, E_nonbond_ww, E_nonbond_qx;
E_nonbonded_t *EQ_nonbond_qq, *EQ_nonbond_qp, *EQ_nonbond_qw, *EQ_nonbond_qx;
E_restraint_t E_restraint, *EQ_restraint;

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

double Ndegf, Ndegfree, Ndegf_solvent, Ndegf_solute, Ndegfree_solvent, Ndegfree_solute;
double Tscale_solute = 1, Tscale_solvent = 1;

void calc_temperature() {
    printf("Ndegf = %f, Ndegfree = %f, n_excluded = %d, Ndegfree_solvent = %f, Ndegfree_solute = %f\n", Ndegf, Ndegfree, n_excluded, Ndegfree_solvent, Ndegfree_solute);
    Temp = 0;
    Tfree = 0;
    double Temp_solute = 0, Tfree_solute = 0, Texcl_solute = 0;
    double Tfree_solvent = 0, Temp_solvent = 0, Texcl_solvent = 0;
    double Ekinmax = 1000.0 * Ndegf * Boltz * md.temperature / 2.0 / n_atoms;
    double ener;
    double mass_i;

    Temp = 0;
    for (int i = 0; i < n_atoms_solute; i++) {
        mass_i = catypes[atypes[i].code - 1].m;
        ener = .5 * mass_i * (pow(velocities[i].x, 2) + pow(velocities[i].y, 2) + pow(velocities[i].z, 2));
        Temp_solute += ener;
        if (!excluded[i]) {
            Tfree_solute += ener;
        }
        else {
            Texcl_solute += ener;
        }
        if (ener > Ekinmax) {
            printf(">>> WARNING: hot atom %d: %f\n", i, ener/Boltz/3);
        }
    }

    for (int i = n_atoms_solute; i < n_atoms; i++) {
        mass_i = catypes[atypes[i].code - 1].m;
        ener = .5 * mass_i * (pow(velocities[i].x, 2) + pow(velocities[i].y, 2) + pow(velocities[i].z, 2));
        Temp_solvent += ener;
        if (!excluded[i]) {
            Tfree_solvent += ener;
        }
        else {
            Texcl_solvent += ener;
        }
        if (ener > Ekinmax) {
            printf(">>> WARNING: hot atom %d: %f\n", i, ener/Boltz/3);
        }
    }

    Tfree = Tfree_solute + Tfree_solvent;
    Temp = Temp_solute + Temp_solvent;

    E_total.Ukin = Temp;

    Temp = 2.0 * Temp / Boltz / Ndegf;
    Tfree = 2.0 * Tfree / Boltz / Ndegfree;

    if (separate_scaling) {
        if (Tfree_solvent != 0) Tscale_solvent = sqrt(1 + (dt / tau_T) * (md.temperature / Tfree_solvent - 1.0));
        if (Tfree_solute != 0) Tscale_solute = sqrt(1 + (dt / tau_T) * (md.temperature / Tfree_solute - 1.0));
    }
    else {
        if (Tfree != 0) Tscale_solvent = sqrt(1 + (dt / tau_T) * (md.temperature / Tfree - 1.0));
        Tscale_solute = Tscale_solvent;
    }
    printf("Tscale = %f, tau_T = %f, Temp = %f, Tfree = %f\n", Tscale_solvent, tau_T, Temp, Tfree);
}

/* =============================================
 * == INTEGRATION METHODS
 * =============================================
 */

void calc_leapfrog() {
    double mass_i, winv_i;
    for (int i = 0; i < n_atoms_solute; i++) {
        mass_i = catypes[atypes[i].code - 1].m;

        winv_i = 1/mass_i;
        velocities[i].x = (velocities[i].x - dvelocities[i].x * dt * winv_i) * Tscale_solute;
        velocities[i].y = (velocities[i].y - dvelocities[i].y * dt * winv_i) * Tscale_solute;
        velocities[i].z = (velocities[i].z - dvelocities[i].z * dt * winv_i) * Tscale_solute;

        // Prepare copy for shake
        xcoords[i].x = coords[i].x;
        xcoords[i].y = coords[i].y;
        xcoords[i].z = coords[i].z;

        coords[i].x += velocities[i].x * dt;
        coords[i].y += velocities[i].y * dt;
        coords[i].z += velocities[i].z * dt;

    }

    for (int i = n_atoms_solute; i < n_atoms; i++) {
        mass_i = catypes[atypes[i].code - 1].m;

        winv_i = 1/mass_i;
        velocities[i].x = (velocities[i].x - dvelocities[i].x * dt * winv_i) * Tscale_solvent;
        velocities[i].y = (velocities[i].y - dvelocities[i].y * dt * winv_i) * Tscale_solvent;
        velocities[i].z = (velocities[i].z - dvelocities[i].z * dt * winv_i) * Tscale_solvent;

        // Prepare copy for shake
        xcoords[i].x = coords[i].x;
        xcoords[i].y = coords[i].y;
        xcoords[i].z = coords[i].z;

        coords[i].x += velocities[i].x * dt;
        coords[i].y += velocities[i].y * dt;
        coords[i].z += velocities[i].z * dt;

    }


    // Shake if necessary
    if (n_shake_constraints > 0) {
        calc_shake_constraints(coords, xcoords);
        for (int i = 0; i < n_atoms; i++) {
            velocities[i].x = (coords[i].x - xcoords[i].x) / dt;
            velocities[i].y = (coords[i].y - xcoords[i].y) / dt;
            velocities[i].z = (coords[i].z - xcoords[i].z) / dt;
        }
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

// Write header (number of atoms & lambdas) to output file
void write_energy_header() {
    FILE * fp;

    char path[1024];
    sprintf(path, "%s/output/%s", base_folder, "energies.csv");

    fp = fopen(path, "w");
  
    fprintf(fp, "%d\n", n_atoms);

    fprintf(fp, "lambdas\n");
    fprintf(fp, "%d\n", n_lambdas);

    for (int state = 0; state < n_lambdas; state++) {
        fprintf(fp, "%f\n", lambdas[state]);
    }
  
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

    fprintf(fp, "interval %d\n", iteration / md.energy);

    fprintf(fp, "[temperature]\n");
    fprintf(fp, "Temp\t%f\n", Temp);
    fprintf(fp, "\n");

    fprintf(fp, "[bonded]\n");
    fprintf(fp, "type\tUbond\tUangle\tUtor\tUimp\n");
    fprintf(fp, "p\t%f\t%f\t%f\t%f\n", E_bond_p.Ubond, E_bond_p.Uangle, E_bond_p.Utor, E_bond_p.Uimp);
    fprintf(fp, "w\t%f\t%f\t%f\t%f\n", E_bond_w.Ubond, E_bond_w.Uangle, E_bond_w.Utor, E_bond_w.Uimp);
    fprintf(fp, "qp\t%f\t%f\t%f\t%f\n", E_bond_q.Ubond, E_bond_q.Uangle, E_bond_q.Utor, E_bond_q.Uimp);
    fprintf(fp, "\n");

    fprintf(fp, "[nonbonded]\n");
    fprintf(fp, "type\tUcoul\tUvdw\n");
    fprintf(fp, "pp\t%f\t%f\n", E_nonbond_pp.Ucoul, E_nonbond_pp.Uvdw);
    fprintf(fp, "pw\t%f\t%f\n", E_nonbond_pw.Ucoul, E_nonbond_pw.Uvdw);
    fprintf(fp, "ww\t%f\t%f\n", E_nonbond_ww.Ucoul, E_nonbond_ww.Uvdw);
    fprintf(fp, "qx\t%f\t%f\n", E_nonbond_qx.Ucoul, E_nonbond_qx.Uvdw);
    fprintf(fp, "\n");

    fprintf(fp, "[restraint]\n");
    fprintf(fp, "Uradx\t%f\n", E_restraint.Uradx);
    fprintf(fp, "Upolx\t%f\n", E_restraint.Upolx);
    fprintf(fp, "Ushell\t%f\n", E_restraint.Ushell);
    fprintf(fp, "Ufix\t%f\n", E_restraint.Ufix);
    fprintf(fp, "Upres\t%f\n", E_restraint.Upres);
    fprintf(fp, "Total\t%f\n", E_restraint.Urestr);
    fprintf(fp, "\n");

    fprintf(fp, "[q-energies]\n");
    fprintf(fp, "lambda\tSUM\tUbond\tUangle\tUtor\tUimp\tUcoul\tUvdw\tUrestr\n");
    for (int state = 0; state < n_lambdas; state++) {
        fprintf(fp, "%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\n", lambdas[state], EQ_total[state].Utot, EQ_bond[state].Ubond,
            EQ_bond[state].Uangle, EQ_bond[state].Utor, EQ_bond[state].Uimp, EQ_nonbond_qx[state].Ucoul, EQ_nonbond_qx[state].Uvdw, EQ_restraint[state].Urestr);
    }
    fprintf(fp, "\n");
    
    fprintf(fp, "[total]\n");
    fprintf(fp, "Ukin\t%f\n", E_total.Ukin);
    fprintf(fp, "Upot\t%f\n", E_total.Upot);
    fprintf(fp, "Utot\t%f\n", E_total.Utot);
    fprintf(fp, "\n");
  
    fclose (fp);
}

void calc_integration() {
    init_variables();

    for (int i = 0; i <= md.steps; i++) {
        calc_integration_step(i);
    }
    
    clean_variables();
}

void reset_energies() {
    for (int i = 0; i < n_atoms; i++) {
        dvelocities[i].x = 0;
        dvelocities[i].y = 0;
        dvelocities[i].z = 0;
    }
    E_total.Upot = 0;
    E_bond_p.Uangle = 0;
    E_bond_p.Ubond = 0;
    E_bond_p.Utor = 0;
    E_bond_p.Uimp = 0;
    E_bond_w.Uangle = 0;
    E_bond_w.Ubond = 0;
    E_bond_w.Utor = 0;
    E_bond_w.Uimp = 0;
    E_bond_q.Uangle = 0;
    E_bond_q.Ubond = 0;
    E_bond_q.Utor = 0;
    E_bond_q.Uimp = 0;
    E_nonbond_pp.Ucoul = 0;
    E_nonbond_pp.Uvdw = 0;
    E_nonbond_pw.Ucoul = 0;
    E_nonbond_pw.Uvdw = 0;
    E_nonbond_ww.Ucoul = 0;
    E_nonbond_ww.Uvdw = 0;
    E_nonbond_qx.Ucoul = 0;
    E_nonbond_qx.Uvdw = 0;
    E_restraint.Uradx = 0;
    E_restraint.Upolx = 0;
    E_restraint.Ufix = 0;
    E_restraint.Ushell = 0;
    E_restraint.Upres = 0;
    E_restraint.Urestr = 0;
    for (int state = 0; state < n_lambdas; state++) {
        EQ_bond[state].Uangle = 0;
        EQ_bond[state].Ubond = 0;
        EQ_bond[state].Utor = 0;
        EQ_bond[state].Uimp = 0;
        EQ_nonbond_qq[state].Ucoul = 0;
        EQ_nonbond_qq[state].Uvdw = 0;
        EQ_nonbond_qp[state].Ucoul = 0;
        EQ_nonbond_qp[state].Uvdw = 0;
        EQ_nonbond_qw[state].Ucoul = 0;
        EQ_nonbond_qw[state].Uvdw = 0;
        EQ_restraint[state].Urestr = 0;
    }
}

void calc_bonded_forces() {
    E_bond_p.Uangle = calc_angle_forces(0, n_angles_solute);
    E_bond_w.Uangle = calc_angle_forces(n_angles_solute, n_angles);

    E_bond_p.Ubond = calc_bond_forces(0, n_bonds_solute);
    E_bond_w.Ubond = calc_bond_forces(n_bonds_solute, n_bonds);

    E_bond_p.Utor = calc_torsion_forces(0, n_torsions_solute);
    E_bond_w.Utor = calc_torsion_forces(n_torsions_solute, n_torsions);

    E_bond_p.Uimp = calc_improper2_forces(0, n_impropers_solute);
    E_bond_w.Uimp = calc_improper2_forces(n_impropers_solute, n_impropers);
}

void calc_integration_step(int iteration) {
    printf("================================================\n");
    if (iteration > 0) {
        printf("== STEP %d\n", iteration);
    }
    else {
        printf("== INITIAL ENERGIES\n");
    }
    printf("================================================\n");

    // Reset derivatives & energies
    reset_energies();

    // Determine temperature and kinetic energy
    calc_temperature();

    // Determine acceleration
    clock_t start = clock();

    // First solute interactions
    calc_bonded_forces();

    clock_t end_bonded = clock();

    clock_t start_pp, end_pp, start_qp, end_qp;
    if (run_gpu) {
        start_qp = clock();
        calc_nonbonded_qp_forces_host();
        end_qp = clock();
        start_pp = clock();
        calc_nonbonded_pp_forces_host();
        end_pp = clock();
    }
    else {
        start_qp = clock();
        calc_nonbonded_qp_forces();
        end_qp = clock();
        start_pp = clock();
        calc_nonbonded_pp_forces();
        end_pp = clock();
    }

    clock_t start_ww, end_ww, start_pw, end_pw;
    // Now solvent interactions
    if (n_waters > 0) {
        if (run_gpu) {
            start_ww = clock();
            calc_nonbonded_ww_forces_host();
            end_ww = clock();
            start_pw = clock();
            calc_nonbonded_pw_forces_host();
            end_pw = clock();
            calc_nonbonded_qw_forces_host();
        }
        else {
            start_ww = clock();
            calc_nonbonded_ww_forces();
            end_ww = clock();
            start_pw = clock();
            calc_nonbonded_pw_forces();
            end_pw = clock();
            calc_nonbonded_qw_forces();
        }
    }

    clock_t end_nonbonded = clock();

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

    calc_restrpos_forces();
    calc_restrang_forces();
    calc_restrwall_forces();
    
    clock_t end_restraints = clock();

    // Q-Q nonbonded interactions
    clock_t start_qq = clock();
    calc_nonbonded_qq_forces();
    clock_t end_qq = clock();

    // Q-atom bonded interactions: loop over Q-atom states
    for (int state = 0; state < n_lambdas; state++) {
        calc_qangle_forces(state);
        calc_qbond_forces(state);
        calc_qtorsion_forces(state);
    }

    clock_t end_qatoms = clock();

    // Now apply leapfrog integration
    calc_leapfrog();

    // Recalculate temperature and kinetic energy for output
    calc_temperature();

    // Update total potential energies with an average of all states
    for (int state = 0; state < n_lambdas; state++) {
        if (lambdas[state] == 0) {
            EQ_bond[state].Uangle = 0;
            EQ_bond[state].Ubond = 0;
            EQ_bond[state].Utor = 0;
            EQ_bond[state].Uimp = 0;
            EQ_nonbond_qq[state].Ucoul = 0;
            EQ_nonbond_qq[state].Uvdw = 0;
            EQ_nonbond_qp[state].Ucoul = 0;
            EQ_nonbond_qp[state].Uvdw = 0;
            EQ_nonbond_qw[state].Ucoul = 0;
            EQ_nonbond_qw[state].Uvdw = 0;
            EQ_restraint[state].Urestr = 0;
        }

        EQ_nonbond_qx[state].Ucoul = EQ_nonbond_qq[state].Ucoul + EQ_nonbond_qp[state].Ucoul + EQ_nonbond_qw[state].Ucoul;
        EQ_nonbond_qx[state].Uvdw = EQ_nonbond_qq[state].Uvdw + EQ_nonbond_qp[state].Uvdw + EQ_nonbond_qw[state].Uvdw;

        EQ_total[state].Utot = EQ_bond[state].Ubond + EQ_bond[state].Uangle + EQ_bond[state].Utor + EQ_bond[state].Uimp
            + EQ_nonbond_qx[state].Ucoul + EQ_nonbond_qx[state].Uvdw + EQ_restraint[state].Urestr;

        E_bond_q.Ubond += EQ_bond[state].Ubond * lambdas[state];
        E_bond_q.Uangle += EQ_bond[state].Uangle * lambdas[state];
        E_bond_q.Utor += EQ_bond[state].Utor * lambdas[state];
        E_bond_q.Uimp += EQ_bond[state].Uimp * lambdas[state];
        E_nonbond_qx.Ucoul += EQ_nonbond_qx[state].Ucoul * lambdas[state];
        E_nonbond_qx.Uvdw += EQ_nonbond_qx[state].Uvdw * lambdas[state];

        // Update protein restraint energies with an average of all states
        E_restraint.Upres += EQ_restraint[state].Urestr * lambdas[state];
    }

    // Update totals
    E_restraint.Urestr = E_restraint.Uradx + E_restraint.Upolx + E_restraint.Ushell + E_restraint.Ufix + E_restraint.Upres;
    E_total.Upot = E_bond_p.Ubond + E_bond_w.Ubond + E_bond_p.Uangle + E_bond_w.Uangle + E_bond_p.Utor + E_bond_p.Uimp
        + E_nonbond_pp.Ucoul + E_nonbond_pp.Uvdw + E_nonbond_pw.Ucoul + E_nonbond_pw.Uvdw + E_nonbond_ww.Ucoul
        + E_nonbond_ww.Uvdw + E_bond_q.Ubond + E_bond_q.Uangle + E_bond_q.Utor + E_bond_q.Uimp
        + E_nonbond_qx.Ucoul + E_nonbond_qx.Uvdw + E_restraint.Urestr;
    E_total.Utot = E_total.Upot + E_total.Ukin;

    printf("[temperature]\n");
    printf("Temp\t%f\n", Temp);
    printf("\n");

    printf("[bonded]\n");
    printf("type\tUbond\tUangle\tUtor\tUimp\n");
    printf("p\t%f\t%f\t%f\t%f\n", E_bond_p.Ubond, E_bond_p.Uangle, E_bond_p.Utor, E_bond_p.Uimp);
    printf("w\t%f\t%f\t%f\t%f\n", E_bond_w.Ubond, E_bond_w.Uangle, E_bond_w.Utor, E_bond_w.Uimp);
    printf("qp\t%f\t%f\t%f\t%f\n", E_bond_q.Ubond, E_bond_q.Uangle, E_bond_q.Utor, E_bond_q.Uimp);
    printf("\n");

    printf("[nonbonded]\n");
    printf("type\tUcoul\tUvdw\n");
    printf("pp\t%f\t%f\n", E_nonbond_pp.Ucoul, E_nonbond_pp.Uvdw);
    printf("pw\t%f\t%f\n", E_nonbond_pw.Ucoul, E_nonbond_pw.Uvdw);
    printf("ww\t%f\t%f\n", E_nonbond_ww.Ucoul, E_nonbond_ww.Uvdw);
    printf("qx\t%f\t%f\n", E_nonbond_qx.Ucoul, E_nonbond_qx.Uvdw);
    printf("\n");

    printf("[restraint]\n");
    printf("Uradx\t%f\n", E_restraint.Uradx);
    printf("Upolx\t%f\n", E_restraint.Upolx);
    printf("Ushell\t%f\n", E_restraint.Ushell);
    printf("Ufix\t%f\n", E_restraint.Ufix);
    printf("Upres\t%f\n", E_restraint.Upres);
    printf("Total\t%f\n", E_restraint.Urestr);
    printf("\n");

    printf("[q-energies]\n");
    printf("lambda\tSUM\tUbond\tUangle\tUtor\tUimp\tUcoul\tUvdw\tUrestr\n");
    for (int state = 0; state < n_lambdas; state++) {
        printf("%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\n", lambdas[state], EQ_total[state].Utot, EQ_bond[state].Ubond,
            EQ_bond[state].Uangle, EQ_bond[state].Utor, EQ_bond[state].Uimp, EQ_nonbond_qx[state].Ucoul, EQ_nonbond_qx[state].Uvdw, EQ_restraint[state].Urestr);
    }
    printf("\n");
    
    printf("[total]\n");
    printf("Ukin\t%f\n", E_total.Ukin);
    printf("Upot\t%f\n", E_total.Upot);
    printf("Utot\t%f\n", E_total.Utot);
    printf("\n");

    // Append output files
    write_coords(iteration);
    write_velocities(iteration);
    write_energies(iteration);    

    clock_t end_calculation = clock();

    // Profiler info
#ifdef __PROFILING__
    printf("Elapsed time for bonded forces: %f\n", (end_bonded-start) / (double)CLOCKS_PER_SEC );
    printf("Elapsed time for non-bonded forces: %f\n", (end_nonbonded-end_bonded) / (double)CLOCKS_PER_SEC);
    printf("Elapsed time for pp interactions: %f\n", (end_pp-start_pp) / (double)CLOCKS_PER_SEC );
    printf("Elapsed time for qq interaction: %f\n",  (end_qq-start_qq) / (double)CLOCKS_PER_SEC );
    printf("Elapsed time for qp interaction: %f\n",  (end_qp-start_qp) / (double)CLOCKS_PER_SEC );
    if (n_waters > 0) {
        printf("Elapsed time for ww interactions: %f\n", (end_ww-start_ww) / (double)CLOCKS_PER_SEC );
        printf("Elapsed time for pw interactions: %f\n", (end_pw-start_pw) / (double)CLOCKS_PER_SEC );
    }
    printf("---\n");
    printf("Elapsed time for entire time-step: %f\n", (end_calculation-start) / (double)CLOCKS_PER_SEC);
#endif /* __PROFILING__ */

}

void init_variables() {
    // From MD file
    init_md("md.csv");

    dt = time_unit * md.stepsize;
    tau_T = time_unit * md.bath_coupling;

    if (run_gpu && n_lambdas > 2) {
        printf(">>> FATAL: More than 2 states not supported on GPU architecture. Exiting...\n");
        exit(EXIT_FAILURE);
    }

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
    init_molecules("molecules.csv");
    init_impropers("impropers.csv");
    init_torsions("torsions.csv");
    init_LJ_matrix();
    init_ngbrs14("ngbrs14.csv");
    init_ngbrs23("ngbrs23.csv");
    init_ngbrs14_long("ngbrs14long.csv");
    init_ngbrs23_long("ngbrs23long.csv");
    // init_restrseqs();
    init_inv_mass();
    if (md.charge_groups) {
        init_charge_groups("charge_groups.csv");
    }

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

    // First part of shrink topology, this needs to be done first as shake constraints are based on bonds
    exclude_qatom_definitions();
    exclude_all_atoms_excluded_definitions();
    
    // Shake constraints, need to be initialized before last part of shrink_topology
    if (md.charge_groups) {
        init_pshells_with_charge_groups();
    }
    else {
        init_pshells();
    }
    init_shake();

    // Now remove shaken bonds
    exclude_shaken_definitions();

    // Init random seed from MD file
    srand(md.random_seed);

    // From calculation in the integration
    init_patoms();
    init_velocities();
    init_dvelocities();
    init_xcoords();

    // From input file
    init_icoords("i_coords.csv");
    init_ivelocities("i_velocities.csv");    

    // Init waters, boundary restrains
    n_waters = (n_atoms - n_atoms_solute) / 3;
    if (n_waters > 0) {
        init_water_sphere();
        init_wshells();
    }

    // Init energy
    EQ_total = (energy_t*) malloc(n_lambdas * sizeof(energy_t));
    EQ_bond = (E_bonded_t*) malloc(n_lambdas * sizeof(E_bonded_t));
    EQ_nonbond_qq = (E_nonbonded_t*) malloc(n_lambdas * sizeof(E_nonbonded_t));
    EQ_nonbond_qp = (E_nonbonded_t*) malloc(n_lambdas * sizeof(E_nonbonded_t));
    EQ_nonbond_qw = (E_nonbonded_t*) malloc(n_lambdas * sizeof(E_nonbonded_t));
    EQ_nonbond_qx = (E_nonbonded_t*) malloc(n_lambdas * sizeof(E_nonbonded_t));
    EQ_restraint = (E_restraint_t*) malloc(n_lambdas * sizeof(E_restraint_t));

    if (n_shake_constraints > 0) {
        initial_shaking();
        stop_cm_translation();
    }
    
    // Write header to file
    write_header("coords.csv");
    write_header("velocities.csv");
    write_energy_header();
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
    free(LJ_matrix);
    free(molecules);
    for (int i = 0; i < n_cgrps_solute + n_cgrps_solvent; i++) {
        free(charge_groups[i].a);
    }
    free(charge_groups);

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
    free(q_atypes);
    free(q_charges);
    free(q_angles);
    free(q_bonds);
    free(q_elscales);
    free(q_exclpairs);
    free(q_impropers);
    free(q_shakes);
    free(q_softcores);
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

    // Shake
    free(mol_n_shakes);    
    free(shake_bonds);
    if (run_gpu) {
        clean_d_solvent();
        clean_d_qatoms();
        clean_d_patoms();
    }

    // Energies & temperature
    free(EQ_total);
    free(EQ_bond);
    free(EQ_nonbond_qq);
    free(EQ_nonbond_qp);
    free(EQ_nonbond_qw);
    free(EQ_nonbond_qx);
    free(EQ_restraint);
}
