#include "system.h"
#include "utils.h"
#include "parse.h"
#include "bonded.h"
#include "nonbonded.h"
#include "solvent.h"
#include "restraints.h"

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
int n_waters;

char base_folder[1024];

double dt = time_unit * step_size;
double tau_T = time_unit * c_bath_coupling;
//double dt = step_size_unit * step_size;

/* =============================================
 * == FROM MD FILE
 * =============================================
 */

md_t md;

void init_md(char *filename) {
    csvfile_t file = read_csv(filename, 0, base_folder);
    char *eptr;

    // [MD]
    md.steps = atoi(file.buffer[1][0]);
    printf("read %d into steps (%s in file)\n", md.steps, file.buffer[1][1]);

    md.stepsize = strtod(file.buffer[2][0], &eptr);
    printf("read %f into stepsize (%s in file)\n", md.stepsize, file.buffer[2][1]);

    md.temperature = strtod(file.buffer[3][0], &eptr);
    printf("read %f into temperature (%s in file)\n", md.temperature, file.buffer[3][1]);

    strcpy(md.thermostat, file.buffer[4][0]);
    printf("read %s into thermostat (%s in file)\n", md.thermostat, file.buffer[4][1]);
    
    md.random_seed = atoi(file.buffer[5][0]);
    printf("read %d into random_seed (%s in file)\n", md.random_seed, file.buffer[5][1]);

    md.initial_temperature = strtod(file.buffer[6][0], &eptr);
    printf("read %f into initial_temperature (%s in file)\n", md.initial_temperature, file.buffer[6][1]);

    md.shake_solvent = file.buffer[7][0] == "on";
    printf("read %s into shake_solvent (%s in file)\n", file.buffer[7][0], file.buffer[7][1]);

    md.shake_hydrogens = file.buffer[8][0] == "on";
    printf("read %s into shake_hydrogens (%s in file)\n", file.buffer[8][0], file.buffer[8][1]);

    md.lrf = file.buffer[9][0] == "on";
    printf("read %s into shake_hydrogens (%s in file)\n", file.buffer[9][0], file.buffer[9][1]);

    // [cut-offs]
    md.solute_solute = strtod(file.buffer[10][0], &eptr);
    printf("read %f into solute_solute (%s in file)\n", md.solute_solute, file.buffer[10][1]);

    md.solvent_solvent = strtod(file.buffer[11][0], &eptr);
    printf("read %f into solvent_solvent (%s in file)\n", md.solvent_solvent, file.buffer[11][1]);

    md.solute_solvent = strtod(file.buffer[12][0], &eptr);
    printf("read %f into solute_solvent (%s in file)\n", md.solute_solvent, file.buffer[12][1]);

    md.q_atom = strtod(file.buffer[13][0], &eptr);
    printf("read %f into q_atom (%s in file)\n", md.q_atom, file.buffer[13][1]);

    // [sphere]
    md.shell_radius = strtod(file.buffer[14][0], &eptr);
    printf("read %f into shell_radius (%s in file)\n", md.shell_radius, file.buffer[14][1]);

    md.shell_force = strtod(file.buffer[15][0], &eptr);
    printf("read %f into shell_force (%s in file)\n", md.shell_force, file.buffer[15][1]);

    // [solvent]
    md.radial_force = strtod(file.buffer[16][0], &eptr);
    printf("read %f into radial_force (%s in file)\n", md.radial_force, file.buffer[16][1]);

    md.polarization = file.buffer[17][0] == "on";
    printf("read %s into polarization (%s in file)\n", file.buffer[17][0], file.buffer[17][1]);

    // [intervals]
    md.non_bond = atoi(file.buffer[18][0]);
    printf("read %d into non_bond (%s in file)\n", md.non_bond, file.buffer[18][1]);

    md.output = atoi(file.buffer[19][0]);
    printf("read %d into output (%s in file)\n", md.output, file.buffer[19][1]);

    md.energy = atoi(file.buffer[20][0]);
    printf("read %d into energy (%s in file)\n", md.energy, file.buffer[20][1]);

    md.trajectory = atoi(file.buffer[21][0]);
    printf("read %d into trajectory (%s in file)\n", md.trajectory, file.buffer[21][1]);

    // [trajectory_atoms]

    // From here on, need a variable to keep track of index in csvfile
    int k = 22;

    // [lambdas]
    n_lambdas = atoi(file.buffer[k][0]);
    printf("reading in %d lambdas\n", n_lambdas);
    lambdas = (double*) malloc(n_lambdas * sizeof(double));
    k++;
    for (int i = 0; i < n_lambdas; i++) {
        lambdas[i] = strtod(file.buffer[k][0], &eptr);
        k++;
    }

    // [sequence_restraints]
    n_restrseqs = atoi(file.buffer[k][0]);
    printf("reading in %d sequence restraints\n", n_restrseqs);
    restrseqs = (restrseq_t*) malloc(n_restrseqs * sizeof(restrseq_t));
    k++;
    for (int i = 0; i < n_restrseqs; i++) {
        restrseq_t restrseq;

        restrseq.ai = atoi(file.buffer[k][0]);
        restrseq.aj = atoi(file.buffer[k][1]);
        restrseq.k = strtod(file.buffer[k][2], &eptr);
        restrseq.ih = file.buffer[k][3] == "1";
        restrseq.to_center = atoi(file.buffer[k][4]);

        restrseqs[i] = restrseq;
        k++;
    }

    // [position_restraints]
    n_restrspos = atoi(file.buffer[k][0]);
    printf("reading in %d position restraints\n", n_restrspos);
    restrspos = (restrpos_t*) malloc(n_restrspos * sizeof(restrpos_t));
    k++;
    for (int i = 0; i < n_restrspos; i++) {
        restrpos_t restrpos;

        restrpos.a = atoi(file.buffer[k][0]);
        restrpos.ipsi = atoi(file.buffer[k][1]);

        coord_t r_x, r_k;

        r_x.x = strtod(file.buffer[k][2], &eptr);
        r_x.y = strtod(file.buffer[k][3], &eptr);
        r_x.z = strtod(file.buffer[k][4], &eptr);
        r_k.x = strtod(file.buffer[k][5], &eptr);
        r_k.y = strtod(file.buffer[k][6], &eptr);
        r_k.z = strtod(file.buffer[k][7], &eptr);

        restrpos.x = r_x;
        restrpos.k = r_k;
        
        restrspos[i] = restrpos;
        k++;
    }

    // [distance_restraints]
    n_restrdists = atoi(file.buffer[k][0]);
    restrdists = (restrdis_t*) malloc(n_restrdists * sizeof(restrdis_t));
    printf("reading in %d distance restraints\n", n_restrdists);
    k++;
    for (int i = 0; i < n_restrdists; i++) {
        restrdis_t restrdist;

        restrdist.ai = atoi(file.buffer[k][0]);
        restrdist.aj = atoi(file.buffer[k][1]);
        restrdist.d1 = strtod(file.buffer[k][2], &eptr);
        restrdist.d2 = strtod(file.buffer[k][3], &eptr);
        restrdist.k = strtod(file.buffer[k][4], &eptr);
        strcpy(restrdist.itext, file.buffer[k][5]);
        strcpy(restrdist.jtext, file.buffer[k][6]);

        restrdists[i] = restrdist;
        k++;
    }

    // [angle_restraints]
    n_restrangs = atoi(file.buffer[k][0]);
    restrangs = (restrang_t*) malloc(n_restrangs * sizeof(restrang_t));
    printf("reading in %d angle restraints\n", n_restrangs);
    k++;
    for (int i = 0; i < n_restrangs; i++) {
        restrang_t restrang;

        restrang.ai = atoi(file.buffer[k][0]);
        restrang.aj = atoi(file.buffer[k][1]);
        restrang.ak = atoi(file.buffer[k][2]);
        restrang.ipsi = atoi(file.buffer[k][3]);
        restrang.ang = strtod(file.buffer[k][4], &eptr);
        restrang.k = strtod(file.buffer[k][5], &eptr);

        restrangs[i] = restrang;
        k++;
    }

    // [wall_restraints]
    n_restrwalls = atoi(file.buffer[k][0]);
    restrwalls = (restrwall_t*) malloc(n_restrwalls * sizeof(restrwall_t));
    printf("reading in %d wall restraints\n", n_restrwalls);
    k++;
    for (int i = 0; i < n_restrwalls; i++) {
        restrwall_t restrwall;

        restrwall.ai = atoi(file.buffer[k][0]);
        restrwall.aj = atoi(file.buffer[k][1]);
        restrwall.d = atoi(file.buffer[k][2]);
        restrwall.k = atoi(file.buffer[k][3]);
        restrwall.dMorse = atoi(file.buffer[k][4]);
        restrwall.aMorse = atoi(file.buffer[k][5]);
        restrwall.ih = file.buffer[k][6] == "1";
        
        restrwalls[i] = restrwall;
        k++;
    }

    clean_csv(file);
}

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

void init_coords(char* filename) {
    csvfile_t file = read_csv(filename, 1, base_folder);

    n_coords = 0;
    n_atoms = 0;
    n_atoms_solute = 0;

    if (file.n_lines < 1) {
        return;
    }

    n_coords = atoi(file.buffer[0][0]);
    n_atoms = n_coords;

    n_atoms_solute = atoi(file.buffer[1][0]);

    coords = (coord_t*) malloc(n_atoms * sizeof(coord_t));
    coords_top = (coord_t*) malloc(n_atoms * sizeof(coord_t));

    for (int i = 0; i < file.n_lines; i++) {
        coord_t coord, coord_top;
        char *eptr;

        coord.x = strtod(file.buffer[i+2][0], &eptr);
        coord.y = strtod(file.buffer[i+2][1], &eptr);
        coord.z = strtod(file.buffer[i+2][2], &eptr);

        coord_top.x = coord.x;
        coord_top.y = coord.y;
        coord_top.z = coord.z;

        coords[i] = coord;
        coords_top[i] = coord_top;
    }
    
    clean_csv(file);
}

void init_bonds(char* filename) {
    csvfile_t file = read_csv(filename, 1, base_folder);
    
    n_bonds = 0;
    n_bonds_solute = 0;

    if (file.n_lines < 1) {
        return;
    }

    n_bonds = atoi(file.buffer[0][0]);
    n_bonds_solute = atoi(file.buffer[1][0]);

    bonds = (bond_t*) malloc(n_bonds * sizeof(bond_t));

    for (int i = 0; i < n_bonds; i++) {
        bond_t bond;

        bond.ai = atoi(file.buffer[i+2][0]);
        bond.aj = atoi(file.buffer[i+2][1]);
        bond.code = atoi(file.buffer[i+2][2]);

        bonds[i] = bond;
    }

    clean_csv(file);
}

void init_cbonds(char* filename) {
    csvfile_t file = read_csv(filename, 0, base_folder);

    n_cbonds = 0;

    if (file.n_lines < 1) {
        return;
    }

    n_cbonds = atoi(file.buffer[0][0]);
    cbonds = (cbond_t*) malloc(n_cbonds * sizeof(cbond_t));

    for (int i = 0; i < n_cbonds; i++) {
        cbond_t cbond;
        char *eptr;

        cbond.code = atoi(file.buffer[i+1][0]);
        cbond.kb = strtod(file.buffer[i+1][1], &eptr);
        cbond.b0 = strtod(file.buffer[i+1][2], &eptr);

        cbonds[i] = cbond;
    }

    clean_csv(file);
}

void init_angles(char* filename) {
    csvfile_t file = read_csv(filename, 1, base_folder);

    n_angles = 0;
    n_angles_solute = 0;

    if (file.n_lines < 1) {
        return;
    }

    n_angles = atoi(file.buffer[0][0]);
    n_angles_solute = atoi(file.buffer[1][0]);

    angles = (angle_t*) malloc(n_angles * sizeof(angle_t));

    for (int i = 0; i < n_angles; i++) {
        angle_t angle;

        angle.ai = atoi(file.buffer[i+2][0]);
        angle.aj = atoi(file.buffer[i+2][1]);
        angle.ak = atoi(file.buffer[i+2][2]);
        angle.code = atoi(file.buffer[i+2][3]);

        angles[i] = angle;
    }

    clean_csv(file);
}

void init_cangles(char* filename) {
    csvfile_t file = read_csv(filename, 0, base_folder);

    n_cangles = 0;
    
    if (file.n_lines < 1) {
        return;
    }

    n_cangles = atoi(file.buffer[0][0]);
    cangles = (cangle_t*) malloc(n_cangles * sizeof(cangle_t));

    for (int i = 0; i < n_cangles; i++) {
        cangle_t cangle;
        char* eptr;

        cangle.code = atoi(file.buffer[i+1][0]);
        cangle.kth = strtod(file.buffer[i+1][1], &eptr);
        cangle.th0 = strtod(file.buffer[i+1][2], &eptr);

        cangles[i] = cangle;
    }
    
    clean_csv(file);
}

void init_excluded(char *filename) {
    excluded = (bool*) calloc(n_atoms, sizeof(bool));
}

void init_torsions(char* filename) {
    csvfile_t file = read_csv(filename, 1, base_folder);

    n_torsions = 0;
    n_torsions_solute = 0;

    if (file.n_lines < 1) {
        return;
    }

    n_torsions = atoi(file.buffer[0][0]);
    n_torsions_solute = atoi(file.buffer[1][0]);

    torsions = (torsion_t*) malloc(n_torsions * sizeof(torsion_t));

    for (int i = 0; i < n_torsions; i++) {
        torsion_t torsion;

        torsion.ai = atoi(file.buffer[i+2][0]);
        torsion.aj = atoi(file.buffer[i+2][1]);
        torsion.ak = atoi(file.buffer[i+2][2]);
        torsion.al = atoi(file.buffer[i+2][3]);
        torsion.code = atoi(file.buffer[i+2][4]);

        torsions[i] = torsion;
    }
    
    clean_csv(file);
}

void init_ctorsions(char* filename) {
    csvfile_t file = read_csv(filename, 0, base_folder);
    
    n_ctorsions = 0;

    if (file.n_lines < 1) {
        return;
    }
    
    n_ctorsions = atoi(file.buffer[0][0]);

    ctorsions = (ctorsion_t*) malloc(n_ctorsions * sizeof(ctorsion_t));

    for (int i = 0; i < n_ctorsions; i++) {
        ctorsion_t ctorsion;
        char* eptr;

        ctorsion.code = atoi(file.buffer[i+1][0]);
        ctorsion.k = strtod(file.buffer[i+1][1], &eptr);
        ctorsion.n = strtod(file.buffer[i+1][2], &eptr);
        ctorsion.d = strtod(file.buffer[i+1][3], &eptr);

        ctorsions[i] = ctorsion;
    }

    clean_csv(file);
}

void init_impropers(char* filename) {
    csvfile_t file = read_csv(filename, 1, base_folder);

    n_impropers = 0;
    n_impropers_solute = 0;

    if (file.n_lines < 1) {
        return;
    }

    n_impropers = atoi(file.buffer[0][0]);
    n_impropers_solute = atoi(file.buffer[1][0]);

    impropers = (improper_t*) malloc(n_impropers * sizeof(improper_t));

    for (int i = 0; i < n_impropers; i++) {
        improper_t improper;

        improper.ai = atoi(file.buffer[i+2][0]);
        improper.aj = atoi(file.buffer[i+2][1]);
        improper.ak = atoi(file.buffer[i+2][2]);
        improper.al = atoi(file.buffer[i+2][3]);
        improper.code = atoi(file.buffer[i+2][4]);

        impropers[i] = improper;
    }
    
    clean_csv(file);
}

void init_cimpropers(char* filename) {
    csvfile_t file = read_csv(filename, 0, base_folder);

    n_cimpropers = 0;

    if (file.n_lines < 1) {
        return;
    }

    n_cimpropers = atoi(file.buffer[0][0]);

    cimpropers = (cimproper_t*) malloc(n_cimpropers * sizeof(cimproper_t));

    for (int i = 0; i < n_cimpropers; i++) {
        cimproper_t cimproper;
        char* eptr;

        cimproper.code = atoi(file.buffer[i+1][0]);
        cimproper.k = strtod(file.buffer[i+1][1], &eptr);
        cimproper.n = strtod(file.buffer[i+1][2], &eptr);
        cimproper.d = strtod(file.buffer[i+1][3], &eptr);

        cimpropers[i] = cimproper;
    }
    
    clean_csv(file);
}

void init_charges(char* filename) {
    csvfile_t file = read_csv(filename, 0, base_folder);

    n_charges = 0;

    if (file.n_lines < 1) {
        return;
    }

    n_charges = atoi(file.buffer[0][0]);

    charges = (charge_t*) malloc(n_charges * sizeof(charge_t));

    for (int i = 0; i < n_charges; i++) {
        charge_t charge;

        charge.a = atoi(file.buffer[i+1][0]);
        charge.code = atoi(file.buffer[i+1][1]);

        charges[i] = charge;
    }

    clean_csv(file);
}

void init_ccharges(char* filename) {
    csvfile_t file = read_csv(filename, 0, base_folder);

    if (file.n_lines < 1) {
        return;
    }
    
    n_ccharges = atoi(file.buffer[0][0]);

    ccharges = (ccharge_t*) malloc(n_ccharges * sizeof(ccharge_t));
    
    for (int i = 0; i < n_ccharges; i++) {
        ccharge_t ccharge;
        char* eptr;

        ccharge.code = atoi(file.buffer[i+1][0]);
        ccharge.charge = strtod(file.buffer[i+1][1], &eptr);

        ccharges[i] = ccharge;
    }

    clean_csv(file);
}

void init_ngbrs14(char* filename) {
    FILE * fp;

    char path[1024];
    sprintf(path, "%s/%s", base_folder, filename);

    if(access(path, F_OK) == -1) {
        printf(">>> FATAL: The following file could not be found. Exiting...\n");
        puts(path);
        exit(EXIT_FAILURE);
    }

    fp = fopen(path, "r");
    if (fp == NULL)
        exit(EXIT_FAILURE);

    char line[1024];
    
    n_ngbrs14 = 0;
    int lines = 0;

    if (fgets(line, 1024, fp)) {
        lines = atoi(line);
    }
    else {
        return;
    }

    ngbrs14 = (ngbr14_t*) malloc(lines * sizeof(ngbr14_t) * line_width);
    int lineI = 0;

    while (fgets(line, 1024, fp)) {
        for (int i = 0; i < line_width; i++) {
            if (line[i] == '1') {
                ngbr14_t ngbr14;
                ngbr14.ai = lineI + 1;
                ngbr14.aj = ((lineI + i + 1) % lines) + 1;
                ngbrs14[n_ngbrs14] = ngbr14;
                n_ngbrs14++;
            }
        }
        lineI++;
    }

    fclose(fp);
}

void init_ngbrs23(char* filename) {
    FILE * fp;

    char path[1024];
    sprintf(path, "%s/%s", base_folder, filename);

    if(access(path, F_OK) == -1) {
        printf(">>> FATAL: The following file could not be found. Exiting...\n");
        puts(path);
        exit(EXIT_FAILURE);
    }

    fp = fopen(path, "r");
    if (fp == NULL)
        exit(EXIT_FAILURE);

    char line[1024];
    
    n_ngbrs23 = 0;
    int lines = 0;

    if (fgets(line, 1024, fp)) {
        lines = atoi(line);
    }
    else {
        return;
    }

    ngbrs23 = (ngbr23_t*) malloc(lines * sizeof(ngbr23_t) * line_width);
    int lineI = 0;

    while (fgets(line, 1024, fp)) {
        for (int i = 0; i < line_width; i++) {
            if (line[i] == '1') {
                ngbr23_t ngbr23;
                ngbr23.ai = lineI + 1;
                ngbr23.aj = ((lineI + i + 1) % lines) + 1;
                ngbrs23[n_ngbrs23] = ngbr23;
                n_ngbrs23++;
            }
        }
        lineI++;
    }

    fclose(fp);
}

void init_catypes(char* filename) {
    csvfile_t file = read_csv(filename, 0, base_folder);

    int n_catomtypes = 0;

    if (file.n_lines < 1) {
        return;
    }

    n_catomtypes = atoi(file.buffer[0][0]);
    catypes = (catype_t*) malloc(n_catomtypes * sizeof(catype_t));

    for (int i = 0; i < n_catomtypes; i++) {
        catype_t catype;
        char* eptr;

        catype.code = atoi(file.buffer[i+1][0]);
        catype.m = strtod(file.buffer[i+1][1], &eptr);
        catype.aii_normal = strtod(file.buffer[i+1][2], &eptr);
        catype.bii_normal = strtod(file.buffer[i+1][3], &eptr);
        catype.aii_polar = strtod(file.buffer[i+1][4], &eptr);
        catype.bii_polar = strtod(file.buffer[i+1][5], &eptr);
        catype.aii_1_4 = strtod(file.buffer[i+1][6], &eptr);
        catype.bii_1_4 = strtod(file.buffer[i+1][7], &eptr);

        catypes[i] = catype;
    }

    clean_csv(file);
}

void init_atypes(char* filename) {
    csvfile_t file = read_csv(filename, 0, base_folder);

    n_atypes = 0;

    if (file.n_lines < 1) {
        return;
    }

    n_atypes = atoi(file.buffer[0][0]);

    atypes = (atype_t*) malloc(n_atypes * sizeof(atype_t));
    for (int i = 0; i < n_atypes; i++) {
        atype_t atype;

        atype.a = atoi(file.buffer[i+1][0]);
        atype.code = atoi(file.buffer[i+1][1]);

        atypes[i] = atype;
    }

    clean_csv(file);
}

/* =============================================
 * == Q ATOMS
 * =============================================
 */

int n_lambdas;
double *lambdas;

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
double A_OO = 0, B_OO, crg_ow, crg_hw, mu_w = 0;

void init_velocities() {
    velocities = (vel_t*) malloc(n_atoms * sizeof(vel_t));

    // If not previous value set, use a Maxwell distribution to fill velocities
    double kT = Boltz * Tmaxw;
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
    Dwmz = 0.26 * exp(-0.19 * (rwater - 15)) + 0.74;
    awmz = 0.2 / (1 + exp(0.4 * (rwater - 25))) + 0.3;

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

    router = rwater;
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
    rin2 = pow(shell_default * rexcl_o, 2);

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
            r2 = pow(coords_top[i].x - centerX, 2) 
                + pow(coords_top[i].y - centerY, 2)
                + pow(coords_top[i].z - centerZ, 2);
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


void init_restrseqs(char* filename) {
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
 * == ENERGY & TEMPERATURE
 * =============================================
 */

void calc_temperature() {
    double Ndegf = 3 * n_atoms;
    double Ekinmax = 1000.0 * Ndegf * Boltz * Temp0 / 2.0 / n_atoms;
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
    
    Tscale = sqrt(1 + (dt / tau_T) * (Temp0 / Temp - 1.0));
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

// Write header (number of atoms) to coordinate output file
void write_header() {
    FILE * fp;

    char path[1024];
    sprintf(path, "%s/output/%s", base_folder, "coords.csv");

    fp = fopen(path, "w");
  
    fprintf(fp, "%d\n", n_atoms);
  
    fclose (fp);
}

// Write step number, coordinates of atoms to coordinate output file
void write_coords(int iteration) {
    if (iteration % 25 != 0) return;
    FILE * fp;
    int i;

    char path[1024];
    sprintf(path, "%s/output/%s", base_folder, "coords.csv");

    fp = fopen(path, "a");

    fprintf(fp, "%d\n", iteration / 25);
    for(i = 0; i < n_atoms; i++) {
        fprintf(fp, "%f;%f;%f\n", coords[i].x, coords[i].y, coords[i].z);
    }
  
    fclose (fp);
}

void calc_integration(int iteration) {
    printf("================================================\n");
    printf("== STEP %d\n", iteration);
    printf("================================================\n");

    // Reset derivatives
    for (int i = 0; i < n_atoms; i++) {
        dvelocities[i].x = 0;
        dvelocities[i].y = 0;
        dvelocities[i].z = 0;
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

    calc_nonbonded_forces();

    clock_t end_nonbonded = clock();

    // Now solvent interactions
    if (n_waters > 0) {
        calc_nonbonded_ww_forces();
        calc_nonbonded_pw_forces();
    }

    // Calculate restraints
    if (n_waters > 0) {
        calc_radix_w_forces();
        calc_polx_w_forces(iteration);
    }
    calc_pshell_forces();
    calc_restrseq_forces();

    // Now apply leapfrog integration
    calc_leapfrog();

    // Recalculate temperature and kinetic energy for output
    calc_temperature();

    // Update totals
    energies.Urestr = energies.Uradx + energies.Upolx + energies.Ushell + energies.Ufix + energies.Upres;
    energies.Upot = energies.Uangle + energies.Ubond + energies.Utor + energies.Ucoul + energies.Uvdw + energies.Urestr;
    energies.Utot = energies.Upot + energies.Ukin;  
    
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

    // Write coordinates to file
    write_coords(iteration);
}


void init_variables() {
    // From topology file
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
    init_restrseqs("restrseqs.csv");

    // From calculation in the integration
    init_velocities();
    init_dvelocities();
    init_xcoords();

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

    // Write header to file
    write_header();
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
    free(ctorsions);
    free(impropers);
    free(torsions);
    free(ngbrs14);
    free(ngbrs23);

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

    // From calculation in the integration
    free(velocities);
    free(dvelocities);
    free(xcoords);
}