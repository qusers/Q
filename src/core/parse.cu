#include "parse.h"
#include "system.h"
#include <stdio.h>
#include <unistd.h>

csvfile_t read_csv(const char* filename, int ext, char* base_folder) {
    csvfile_t retval;

    retval.ext = ext;

    char path[1024];
    sprintf(path, "%s/%s", base_folder, filename);
    if(access(path, F_OK) == -1) {
        printf(">>> FATAL: The following file could not be found. Exiting...\n");
        puts(path);
        exit(EXIT_FAILURE);
    }

    // File handle
    FILE * fp;

    fp = fopen(path, "r");
    if (fp == NULL)
        exit(EXIT_FAILURE);

    // Get number of lines
    char nlines[COLUMN_WIDTH];
    
    if (fgets(nlines, COLUMN_WIDTH, fp)) {
        retval.n_lines = atoi(nlines);
    }
    else {
        retval.n_lines = 0;
        return retval;
    }

    if (retval.n_lines == 0) {
        return retval;
    }

    char line[N_COLUMNS * COLUMN_WIDTH];
    retval.buffer = (char***) malloc(retval.n_lines * N_COLUMNS * sizeof(char**));

    for (int i = 0; i <= retval.n_lines + ext; i++) {
        retval.buffer[i] = (char**) malloc(N_COLUMNS * sizeof(char*));
        for (int j = 0; j < N_COLUMNS; j++) {
            retval.buffer[i][j] = (char*) malloc(COLUMN_WIDTH * sizeof(char));
        }
    }

    strcpy(retval.buffer[0][0], nlines);
    int lineI = 1;

    // Read in file
    while (fgets(line, N_COLUMNS * COLUMN_WIDTH, fp)) {
        int field = 0;
        // NOTE strtok clobbers tmp
        char* tmp = strdup(line);
        const char* tok;
        for (tok = strtok(tmp, ";");
            tok && *tok;
            tok = strtok(NULL, ";\n"))
        {
            strcpy(retval.buffer[lineI][field], tok);
            field++;
        }
        free(tmp);
        lineI++;
    }

    fclose(fp);

    return retval;
}

void clean_csv(csvfile_t file) {
    if (file.n_lines > 0) {
        for (int i = 0; i <= file.n_lines + file.ext; i++) {
            for (int j = 0; j < N_COLUMNS; j++) {
                free(file.buffer[i][j]);
            }
            free(file.buffer[i]);
        }
        free(file.buffer);
    }
}


/* =============================================
 * == FROM MD FILE
 * =============================================
 */

void init_md(const char *filename) {
    csvfile_t file = read_csv(filename, 0, base_folder);
    char *eptr;

    // [MD]
    md.steps = atoi(file.buffer[1][1]);
    #ifdef VERBOSE
    printf("read %d into steps (%s in file)\n", md.steps, file.buffer[1][0]);
    #endif
    md.stepsize = strtod(file.buffer[2][1], &eptr);
    #ifdef VERBOSE
    printf("read %f into stepsize (%s in file)\n", md.stepsize, file.buffer[2][0]);
    #endif
    md.temperature = strtod(file.buffer[3][1], &eptr);
    #ifdef VERBOSE
    printf("read %f into temperature (%s in file)\n", md.temperature, file.buffer[3][0]);
    #endif
    strcpy(md.thermostat, file.buffer[4][1]);
    #ifdef VERBOSE
    printf("read %s into thermostat (%s in file)\n", md.thermostat, file.buffer[4][0]);
    #endif
    md.bath_coupling = strtod(file.buffer[5][1], &eptr);
    #ifdef VERBOSE
    printf("read %f into bath_coupling (%s in file)\n", md.bath_coupling, file.buffer[5][0]);
    #endif
    md.random_seed = atoi(file.buffer[6][1]);
    #ifdef VERBOSE
    printf("read %d into random_seed (%s in file)\n", md.random_seed, file.buffer[6][0]);
    #endif
    md.initial_temperature = strtod(file.buffer[7][1], &eptr);
    #ifdef VERBOSE
    printf("read %f into initial_temperature (%s in file)\n", md.initial_temperature, file.buffer[7][0]);
    #endif
    md.shake_solvent = strcmp(file.buffer[8][1], "on") == 0;
    #ifdef VERBOSE
    printf("read %s into shake_solvent (%s in file)\n", file.buffer[8][1], file.buffer[8][0]);
    #endif
    md.shake_solute = strcmp(file.buffer[9][1], "on") == 0;
    #ifdef VERBOSE
    printf("read %s into shake_solute (%s in file)\n", file.buffer[9][1], file.buffer[9][0]);
    #endif
    md.shake_hydrogens = strcmp(file.buffer[10][1], "on") == 0;
    #ifdef VERBOSE
    printf("read %s into shake_hydrogens (%s in file)\n", file.buffer[10][1], file.buffer[10][0]);
    #endif
    md.lrf = strcmp(file.buffer[11][1], "on") == 0;
    #ifdef VERBOSE
    printf("read %s into lrf (%s in file)\n", file.buffer[11][1], file.buffer[11][0]);
    #endif
    md.charge_groups = strcmp(file.buffer[12][1], "on") == 0;
    #ifdef VERBOSE
    printf("read %s into charge_groups (%s in file)\n", file.buffer[12][1], file.buffer[12][0]);
    #endif
    // [cut-offs]
    md.solute_solute = strtod(file.buffer[13][1], &eptr);
    #ifdef VERBOSE
    printf("read %f into solute_solute (%s in file)\n", md.solute_solute, file.buffer[13][0]);
    #endif
    md.solvent_solvent = strtod(file.buffer[14][1], &eptr);
    #ifdef VERBOSE
    printf("read %f into solvent_solvent (%s in file)\n", md.solvent_solvent, file.buffer[14][0]);
    #endif
    md.solute_solvent = strtod(file.buffer[15][1], &eptr);
    #ifdef VERBOSE
    printf("read %f into solute_solvent (%s in file)\n", md.solute_solvent, file.buffer[15][0]);
    #endif
    md.q_atom = strtod(file.buffer[16][1], &eptr);
    #ifdef VERBOSE
    printf("read %f into q_atom (%s in file)\n", md.q_atom, file.buffer[16][0]);
    #endif
    // [sphere]
    md.shell_radius = strtod(file.buffer[17][1], &eptr);
    #ifdef VERBOSE
    printf("read %f into shell_radius (%s in file)\n", md.shell_radius, file.buffer[17][0]);
    #endif
    md.shell_force = strtod(file.buffer[18][1], &eptr);
    #ifdef VERBOSE
    printf("read %f into shell_force (%s in file)\n", md.shell_force, file.buffer[18][0]);
    #endif
    // [solvent]
    md.radial_force = strtod(file.buffer[19][1], &eptr);
    #ifdef VERBOSE
    printf("read %f into radial_force (%s in file)\n", md.radial_force, file.buffer[19][0]);
    #endif
    md.polarisation = true;
    #ifdef VERBOSE
    printf("read %s into polarisation (%s in file)\n", file.buffer[20][1], file.buffer[20][0]);
    #endif
    md.polarisation_force = strtod(file.buffer[21][1], &eptr);
    #ifdef VERBOSE
    printf("read %s into polarisation_force (%s in file)\n", file.buffer[21][1], file.buffer[21][0]);
    #endif
    // [intervals]
    md.non_bond = atoi(file.buffer[22][1]);
    #ifdef VERBOSE
    printf("read %d into non_bond (%s in file)\n", md.non_bond, file.buffer[22][0]);
    #endif
    md.output = atoi(file.buffer[23][1]);
    #ifdef VERBOSE
    printf("read %d into output (%s in file)\n", md.output, file.buffer[23][0]);
    #endif
    md.energy = atoi(file.buffer[24][1]);
    #ifdef VERBOSE
    printf("read %d into energy (%s in file)\n", md.energy, file.buffer[24][0]);
    #endif
    md.trajectory = atoi(file.buffer[25][1]);
    #ifdef VERBOSE
    printf("read %d into trajectory (%s in file)\n", md.trajectory, file.buffer[25][0]);
    #endif
    // [trajectory_atoms]

    // From here on, need a variable to keep track of index in csvfile
    int k = 26;

    // [lambdas]
    n_lambdas = atoi(file.buffer[k][0]);
    #ifdef VERBOSE
    printf("reading in %d lambdas (%s in file)\n", n_lambdas, file.buffer[k][1]);
    #endif
    lambdas = (double*) malloc(n_lambdas * sizeof(double));
    k++;
    for (int i = 0; i < n_lambdas; i++) {
        lambdas[i] = strtod(file.buffer[k][0], &eptr);
        k++;
    }

    // [sequence_restraints]
    printf("k = %d\n", k);
    n_restrseqs = atoi(file.buffer[k][0]);
    printf("reading in %d sequence restraints (%s in file)\n", n_restrseqs, file.buffer[k][1]);
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
    printf("reading in %d position restraints\n (%s in file )", n_restrspos, file.buffer[k][1]);
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
    printf("reading in %d distance restraints (%s in file)\n", n_restrdists, file.buffer[k][1]);
    k++;
    for (int i = 0; i < n_restrdists; i++) {
        restrdis_t restrdist;

        restrdist.ai = atoi(file.buffer[k][0]);
        restrdist.aj = atoi(file.buffer[k][1]);
        restrdist.d1 = strtod(file.buffer[k][2], &eptr);
        restrdist.d2 = strtod(file.buffer[k][3], &eptr);
        restrdist.k = strtod(file.buffer[k][4], &eptr);
        restrdist.ipsi = atoi(file.buffer[k][5]);

        restrdists[i] = restrdist;
        k++;
    }

    // [angle_restraints]
    n_restrangs = atoi(file.buffer[k][0]);
    restrangs = (restrang_t*) malloc(n_restrangs * sizeof(restrang_t));
    printf("reading in %d angle restraints (%s in file)\n", n_restrangs, file.buffer[k][1]);
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
    printf("reading in %d wall restraints (%s in file)\n", n_restrwalls, file.buffer[k][1]);
    k++;
    for (int i = 0; i < n_restrwalls; i++) {
        restrwall_t restrwall;

        restrwall.ai = atoi(file.buffer[k][0]);
        restrwall.aj = atoi(file.buffer[k][1]);
        restrwall.d = strtod(file.buffer[k][2], &eptr);
        restrwall.k = strtod(file.buffer[k][3], &eptr);
        restrwall.dMorse = strtod(file.buffer[k][4], &eptr);
        restrwall.aMorse = strtod(file.buffer[k][5], &eptr);
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

void init_topo(const char *filename) {
    csvfile_t file = read_csv(filename, 0, base_folder);
    char *eptr;

    coord_t solute_center, solvent_center;

    topo.solvent_type = atoi(file.buffer[1][0]);
    topo.exclusion_radius = strtod(file.buffer[2][0], &eptr);
    topo.solvent_radius = strtod(file.buffer[3][0], &eptr);
    solute_center.x = strtod(file.buffer[4][0], &eptr);
    solute_center.y = strtod(file.buffer[4][1], &eptr);
    solute_center.z = strtod(file.buffer[4][2], &eptr);
    solvent_center.x = strtod(file.buffer[5][0], &eptr);
    solvent_center.y = strtod(file.buffer[5][1], &eptr);
    solvent_center.z = strtod(file.buffer[5][2], &eptr);

    topo.solute_center = solute_center;
    topo.solvent_center = solvent_center;

    topo.el14_scale = strtod(file.buffer[6][0], &eptr);
    topo.coulomb_constant = strtod(file.buffer[7][0], &eptr);

    clean_csv(file);
}

void init_coords(const char* filename) {
    csvfile_t file = read_csv(filename, 1, base_folder);

    n_coords = 0;
    n_atoms = 0;
    n_atoms_solute = 0;

    if (file.n_lines < 1) {
        clean_csv(file);
        return;
    }

    n_coords = atoi(file.buffer[0][0]);
    n_atoms = n_coords;

    n_atoms_solute = atoi(file.buffer[1][0]);

    coords = (coord_t*) malloc(n_atoms * sizeof(coord_t));
    coords_top = (coord_t*) malloc(n_atoms * sizeof(coord_t));

    for (int i = 0; i < file.n_lines; i++) {
        char *eptr;

        coords[i].x = strtod(file.buffer[i+2][0], &eptr);
        coords[i].y = strtod(file.buffer[i+2][1], &eptr);
        coords[i].z = strtod(file.buffer[i+2][2], &eptr);

        coords_top[i].x = coords[i].x;
        coords_top[i].y = coords[i].y;
        coords_top[i].z = coords[i].z;
    }
    
    clean_csv(file);
}

void init_bonds(const char* filename) {
    csvfile_t file = read_csv(filename, 1, base_folder);
    
    n_bonds = 0;
    n_bonds_solute = 0;

    if (file.n_lines < 1) {
        clean_csv(file);
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

void init_cbonds(const char* filename) {
    csvfile_t file = read_csv(filename, 0, base_folder);

    n_cbonds = 0;

    if (file.n_lines < 1) {
        clean_csv(file);
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

void init_angles(const char* filename) {
    csvfile_t file = read_csv(filename, 1, base_folder);

    n_angles = 0;
    n_angles_solute = 0;

    if (file.n_lines < 1) {
        clean_csv(file);
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

void init_cangles(const char* filename) {
    csvfile_t file = read_csv(filename, 0, base_folder);

    n_cangles = 0;
    
    if (file.n_lines < 1) {
        clean_csv(file);
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

void init_excluded(const char *filename) {
    excluded = (bool*) malloc(n_atoms * sizeof(bool));
    n_excluded = 0;

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

    char line[8192];
    
    n_excluded = 0;

    if (fgets(line, 8192, fp)) {
        for (int i = 0; i < n_atoms; i++) {
            bool excl = (line[i] == '1');
            excluded[i] = excl;
            if (excl) {
                n_excluded++;
            }
        }
    }

    fclose(fp);
}

void init_torsions(const char* filename) {
    csvfile_t file = read_csv(filename, 1, base_folder);

    n_torsions = 0;
    n_torsions_solute = 0;

    if (file.n_lines < 1) {
        clean_csv(file);
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

void init_ctorsions(const char* filename) {
    csvfile_t file = read_csv(filename, 0, base_folder);
    
    n_ctorsions = 0;

    if (file.n_lines < 1) {
        clean_csv(file);
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

void init_impropers(const char* filename) {
    csvfile_t file = read_csv(filename, 1, base_folder);

    n_impropers = 0;
    n_impropers_solute = 0;

    if (file.n_lines < 1) {
        clean_csv(file);
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

void init_cimpropers(const char* filename) {
    csvfile_t file = read_csv(filename, 0, base_folder);

    n_cimpropers = 0;

    if (file.n_lines < 1) {
        clean_csv(file);
        return;
    }

    n_cimpropers = atoi(file.buffer[0][0]);

    cimpropers = (cimproper_t*) malloc(n_cimpropers * sizeof(cimproper_t));

    for (int i = 0; i < n_cimpropers; i++) {
        cimproper_t cimproper;
        char* eptr;

        cimproper.code = atoi(file.buffer[i+1][0]);
        cimproper.k = strtod(file.buffer[i+1][1], &eptr);
        cimproper.phi0 = strtod(file.buffer[i+1][2], &eptr);

        cimpropers[i] = cimproper;
    }
    
    clean_csv(file);
}

void init_charges(const char* filename) {
    csvfile_t file = read_csv(filename, 0, base_folder);

    n_charges = 0;

    if (file.n_lines < 1) {
        clean_csv(file);
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

void init_ccharges(const char* filename) {
    csvfile_t file = read_csv(filename, 0, base_folder);

    if (file.n_lines < 1) {
        clean_csv(file);
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

void init_ngbrs14(const char* filename) {
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
    
    int lines = 0;

    if (fgets(line, 1024, fp)) {
        lines = atoi(line);
    }
    else {
        return;
    }

    int lineI = 0;

    while (fgets(line, 1024, fp)) {
        for (int i = 0; i < line_width; i++) {
            if (line[i] == '1') {
                int ix = lineI;
                int jx = (lineI + i + 1) % lines;
                // if (ix < 100 && jx < 100) printf("i = %d j = %d\n", ix+1, jx+1);
                LJ_matrix[ix * n_atoms_solute + jx] = 1;
            }
        }
        lineI++;
    }

    fclose(fp);
}

void init_ngbrs23(const char* filename) {
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
    
    int lines = 0;

    if (fgets(line, 1024, fp)) {
        lines = atoi(line);
    }
    else {
        return;
    }

    int lineI = 0;

    while (fgets(line, 1024, fp)) {
        for (int i = 0; i < line_width; i++) {
            if (line[i] == '1') {
                int ix = lineI;
                int jx = (lineI + i + 1) % lines;
                // if (ix < 100 && jx < 100) printf("i = %d j = %d\n", ix+1, jx+1);
                LJ_matrix[ix * n_atoms_solute + jx] = 3;
            }
        }
        lineI++;
    }

    fclose(fp);
}

void init_ngbrs14_long(const char* filename) {
    csvfile_t file = read_csv(filename, 0, base_folder);

    if (file.n_lines < 1) {
        clean_csv(file);
        return;
    }

    int n_ngbrs14_long = atoi(file.buffer[0][0]);

    for (int i = 0; i < n_ngbrs14_long; i++) {
        int ix = atoi(file.buffer[i+1][0])-1;
        int jx = atoi(file.buffer[i+1][1])-1;
        LJ_matrix[ix * n_atoms_solute + jx] = 1;
    }

    clean_csv(file);
}

void init_ngbrs23_long(const char* filename) {
    csvfile_t file = read_csv(filename, 0, base_folder);

    if (file.n_lines < 1) {
        clean_csv(file);
        return;
    }

    int n_ngbrs23_long = atoi(file.buffer[0][0]);

    for (int i = 0; i < n_ngbrs23_long; i++) {
        int ix = atoi(file.buffer[i+1][0])-1;
        int jx = atoi(file.buffer[i+1][1])-1;
        LJ_matrix[ix * n_atoms_solute + jx] = 3;
    }

    clean_csv(file);
}

void init_LJ_matrix() {
    LJ_matrix = (int *) malloc(n_atoms_solute * n_atoms_solute * sizeof(int));
    memset(LJ_matrix, 0, n_atoms_solute * n_atoms_solute * sizeof(int));
}

void init_catypes(const char* filename) {
    csvfile_t file = read_csv(filename, 0, base_folder);

    if (file.n_lines < 1) {
        clean_csv(file);
        return;
    }

    n_catypes = atoi(file.buffer[0][0]);
    catypes = (catype_t*) malloc(n_catypes * sizeof(catype_t));

    for (int i = 0; i < n_catypes; i++) {
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

void init_atypes(const char* filename) {
    csvfile_t file = read_csv(filename, 0, base_folder);

    n_atypes = 0;

    if (file.n_lines < 1) {
        clean_csv(file);
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

void init_molecules(const char *filename) {
    csvfile_t file = read_csv(filename, 0, base_folder);
    
    if (file.n_lines < 1) {
        return;
    }

    n_molecules = atoi(file.buffer[0][0]);
    molecules = (int*) malloc(n_molecules * sizeof(int));

    for (int i = 0; i < n_molecules; i++) {
        molecules[i] = atoi(file.buffer[i+1][0]);
    }

    clean_csv(file);
}

void init_charge_groups(const char *filename) {
    csvfile_t file = read_csv(filename, 0, base_folder);

    if (file.n_lines < 1) {
        return;
    }

    n_cgrps_solute = atoi(file.buffer[1][0]);
    n_cgrps_solvent = atoi(file.buffer[1][1]);

    int n_charge_groups = n_cgrps_solute + n_cgrps_solvent;

    charge_groups = (cgrp_t*) malloc(n_charge_groups * sizeof(cgrp_t));

    int line_nr = 2;
    int n_atoms_crgp = 0;

    for (int i = 0; i < n_charge_groups; i++) {
        cgrp_t charge_group;

        n_atoms_crgp = atoi(file.buffer[line_nr][0]);
        charge_group.n_atoms = n_atoms_crgp;
        charge_group.iswitch = atoi(file.buffer[line_nr][1]);
        charge_group.a = (int*) malloc(n_atoms_crgp * sizeof(int));
        
        line_nr++;
        for (int j = 0; j < charge_group.n_atoms; j++) {
            charge_group.a[j] = atoi(file.buffer[line_nr][0]);
            line_nr++;
        }

        charge_groups[i] = charge_group;
    }

    clean_csv(file);
}

/* =============================================
 * == FROM FEP FILE
 * =============================================
 */

void init_qangcouples(const char *filename) {
    csvfile_t file = read_csv(filename, 0, base_folder);

    n_qangcouples = 0;

    if (file.n_lines < 1) {
        clean_csv(file);
        return;
    }

    n_qangcouples = atoi(file.buffer[0][0]);
    q_angcouples = (q_angcouple_t*) malloc (n_qangcouples * sizeof(q_angcouple_t));

    for (int i = 0; i < n_qangcouples; i++) {
        q_angcouples[i].acode = atoi(file.buffer[i+1][0]);
        q_angcouples[i].bcode = atoi(file.buffer[i+1][1]);
    }

    clean_csv(file);
}

void init_qatoms(const char *filename) {
    csvfile_t file = read_csv(filename, 0, base_folder);

    n_qatoms = 0;

    if (file.n_lines < 1) {
        clean_csv(file);
        return;
    }

    n_qatoms = atoi(file.buffer[0][0]);
    q_atoms = (q_atom_t*) malloc (n_qatoms * sizeof(q_atom_t));

    for (int i = 0; i < n_qatoms; i++) {
        q_atoms[i].a = atoi(file.buffer[i+1][0]);
    }

    clean_csv(file);
}

void init_qcangles(const char *filename) {
    csvfile_t file = read_csv(filename, 0, base_folder);

    n_qcangles = 0;

    if (file.n_lines < 1) {
        clean_csv(file);
        return;
    }

    n_qcangles = atoi(file.buffer[0][0]);
    q_cangles = (q_cangle_t*) malloc (n_qcangles * sizeof(q_cangle_t));

    for (int i = 0; i < n_qcangles; i++) {
        char *eptr;
        q_cangles[i].kth = strtod(file.buffer[i+1][0], &eptr);
        q_cangles[i].th0 = strtod(file.buffer[i+1][1], &eptr);
    }

    clean_csv(file);
}

void init_qcatypes(const char*filename) {
    csvfile_t file = read_csv(filename, 0, base_folder);

    n_qcatypes = 0;

    if (file.n_lines < 1) {
        clean_csv(file);
        return;
    }

    n_qcatypes = atoi(file.buffer[0][0]);
    q_catypes = (q_catype_t*) malloc (n_qcatypes * sizeof(q_catype_t));

    for (int i = 0; i < n_qcatypes; i++) {
        char *eptr;
        strcpy(q_catypes[i].name, file.buffer[i+1][0]);
        q_catypes[i].Ai = strtod(file.buffer[i+1][1], &eptr);
        q_catypes[i].Bi = strtod(file.buffer[i+1][2], &eptr);
        q_catypes[i].Ci = strtod(file.buffer[i+1][3], &eptr);
        q_catypes[i].ai = strtod(file.buffer[i+1][4], &eptr);
        q_catypes[i].Ai_14 = strtod(file.buffer[i+1][5], &eptr);
        q_catypes[i].Bi_14 = strtod(file.buffer[i+1][6], &eptr);
        q_catypes[i].m = strtod(file.buffer[i+1][7], &eptr);
    }

    clean_csv(file);
}

void init_qcbonds(const char*filename) {
    csvfile_t file = read_csv(filename, 0, base_folder);

    n_qcbonds = 0;

    if (file.n_lines < 1) {
        clean_csv(file);
        return;
    }

    n_qcbonds = atoi(file.buffer[0][0]);
    q_cbonds = (q_cbond_t*) malloc (n_qcbonds * sizeof(q_cbond_t));

    for (int i = 0; i < n_qcbonds; i++) {
        char *eptr;
        q_cbonds[i].kb = strtod(file.buffer[i+1][0], &eptr);
        q_cbonds[i].b0 = strtod(file.buffer[i+1][1], &eptr);
    }

    clean_csv(file);
}

void init_qcimpropers(const char*filename) {
    csvfile_t file = read_csv(filename, 0, base_folder);

    n_qcimpropers = 0;

    if (file.n_lines < 1) {
        clean_csv(file);
        return;
    }

    n_qcimpropers = atoi(file.buffer[0][0]);
    q_cimpropers = (q_cimproper_t*) malloc (n_qcimpropers * sizeof(q_cimproper_t));

    for (int i = 0; i < n_qcimpropers; i++) {
        char *eptr;
        q_cimpropers[i].k = strtod(file.buffer[i+1][0], &eptr);
        q_cimpropers[i].phi0 = strtod(file.buffer[i+1][1], &eptr);
    }

    clean_csv(file);
}

void init_qctorsions(const char*filename) {
    csvfile_t file = read_csv(filename, 0, base_folder);

    n_qctorsions = 0;

    if (file.n_lines < 1) {
        clean_csv(file);
        return;
    }

    n_qctorsions = atoi(file.buffer[0][0]);
    q_ctorsions = (q_ctorsion_t*) malloc (n_qctorsions * sizeof(q_ctorsion_t));

    for (int i = 0; i < n_qctorsions; i++) {
        char *eptr;
        q_ctorsions[i].k = strtod(file.buffer[i+1][0], &eptr);
        q_ctorsions[i].n = strtod(file.buffer[i+1][1], &eptr);
        q_ctorsions[i].d = strtod(file.buffer[i+1][2], &eptr);
    }

    clean_csv(file);
}

void init_qoffdiags(const char*filename) {
    csvfile_t file = read_csv(filename, 0, base_folder);

    n_qoffdiags = 0;

    if (file.n_lines < 1) {
        clean_csv(file);
        return;
    }

    n_qoffdiags = atoi(file.buffer[0][0]);
    q_offdiags = (q_offdiag_t*) malloc (n_qoffdiags * sizeof(q_offdiag_t));

    for (int i = 0; i < n_qoffdiags; i++) {
        char *eptr;
        q_offdiags[i].i = atoi(file.buffer[i+1][0]);
        q_offdiags[i].j = atoi(file.buffer[i+1][1]);
        q_offdiags[i].qk = atoi(file.buffer[i+1][2]);
        q_offdiags[i].ql = atoi(file.buffer[i+1][3]);
        q_offdiags[i].Aij = strtod(file.buffer[i+1][4], &eptr);
        q_offdiags[i].muij = strtod(file.buffer[i+1][5], &eptr);
    }

    clean_csv(file);
}

void init_qimprcouples(const char*filename) {
    csvfile_t file = read_csv(filename, 0, base_folder);

    n_qimprcouples = 0;

    if (file.n_lines < 1) {
        clean_csv(file);
        return;
    }

    n_qimprcouples = atoi(file.buffer[0][0]);
    q_imprcouples = (q_imprcouple_t*) malloc (n_qimprcouples * sizeof(q_imprcouple_t));

    for (int i = 0; i < n_qimprcouples; i++) {
        q_imprcouples[i].icode = atoi(file.buffer[i+1][0]);
        q_imprcouples[i].bcode = atoi(file.buffer[i+1][1]);
    }

    clean_csv(file);
}

void init_qsoftpairs(const char*filename) {
    csvfile_t file = read_csv(filename, 0, base_folder);

    n_qsoftpairs = 0;

    if (file.n_lines < 1) {
        clean_csv(file);
        return;
    }

    n_qsoftpairs = atoi(file.buffer[0][0]);
    q_softpairs = (q_softpair_t*) malloc (n_qsoftpairs * sizeof(q_softpair_t));

    for (int i = 0; i < n_qsoftpairs; i++) {
        q_softpairs[i].qi = atoi(file.buffer[i+1][0]);
        q_softpairs[i].qj = atoi(file.buffer[i+1][1]);
    }

    clean_csv(file);
}

void init_qtorcouples(const char*filename) {
    csvfile_t file = read_csv(filename, 0, base_folder);

    n_qtorcouples = 0;

    if (file.n_lines < 1) {
        clean_csv(file);
        return;
    }

    n_qtorcouples = atoi(file.buffer[0][0]);
    q_torcouples = (q_torcouple_t*) malloc (n_qtorcouples * sizeof(q_torcouple_t));

    for (int i = 0; i < n_qtorcouples; i++) {
        q_torcouples[i].tcode = atoi(file.buffer[i+1][0]);
        q_torcouples[i].bcode = atoi(file.buffer[i+1][1]);
    }

    clean_csv(file);
}

void init_qangles(const char*filename) {
    csvfile_t file = read_csv(filename, 0, base_folder);

    if (file.n_lines < 1) {
        clean_csv(file);
        return;
    }

    n_qangles = atoi(file.buffer[0][0]) / n_lambdas;
    q_angles = (q_angle_t*) malloc (n_qangles * n_lambdas * sizeof(q_angle_t));

    for (int i = 0; i < n_qangles; i++) {
        for (int j = 0; j < n_lambdas; j++) {
            q_angles[i + j * n_qangles].ai = atoi(file.buffer[i + j * n_qangles + 1][0]);
            q_angles[i + j * n_qangles].aj = atoi(file.buffer[i + j * n_qangles + 1][1]);
            q_angles[i + j * n_qangles].ak = atoi(file.buffer[i + j * n_qangles + 1][2]);
            q_angles[i + j * n_qangles].code = atoi(file.buffer[i + j * n_qangles + 1][3]);
        }
    }

    clean_csv(file);
}

void init_qatypes(const char*filename) {
    csvfile_t file = read_csv(filename, 0, base_folder);

    if (file.n_lines < 1) {
        clean_csv(file);
        return;
    }

    q_atypes = (q_atype_t*) malloc (n_qatoms * n_lambdas * sizeof(q_atype_t));
    for (int i = 0; i < n_qatoms; i++) {
        for (int j = 0; j < n_lambdas; j++) {
            q_atypes[i + j * n_qatoms].code = atoi(file.buffer[i + j * n_qatoms + 1][0]);
        }
    }

    clean_csv(file);
}

void init_qbonds(const char*filename) {
    csvfile_t file = read_csv(filename, 0, base_folder);

    if (file.n_lines < 1) {
        clean_csv(file);
        return;
    }

    n_qbonds = atoi(file.buffer[0][0]) / n_lambdas;
    q_bonds = (q_bond_t*) malloc (n_qbonds * n_lambdas * sizeof(q_bond_t));

    for (int i = 0; i < n_qbonds; i++) {
        for (int j = 0; j < n_lambdas; j++) {
            q_bonds[i + j * n_qbonds].ai = atoi(file.buffer[i + j * n_qbonds + 1][0]);
            q_bonds[i + j * n_qbonds].aj = atoi(file.buffer[i + j * n_qbonds + 1][1]);
            q_bonds[i + j * n_qbonds].code = atoi(file.buffer[i + j * n_qbonds + 1][2]);
        }
    }

    clean_csv(file);
}

void init_qcharges(const char*filename) {
    csvfile_t file = read_csv(filename, 0, base_folder);

    if (file.n_lines < 1) {
        clean_csv(file);
        return;
    }

    q_charges = (q_charge_t*) malloc (n_qatoms * n_lambdas * sizeof(q_charge_t));

    for (int i = 0; i < n_qatoms; i++) {
        for (int j = 0; j < n_lambdas; j++) {
            char *eptr;
            q_charges[i + j * n_qatoms].q = strtod(file.buffer[i + j * n_qatoms + 1][0], &eptr);
        }
    }

    clean_csv(file);
}

void init_qelscales(const char*filename) {
    csvfile_t file = read_csv(filename, 0, base_folder);

    if (file.n_lines < 1) {
        clean_csv(file);
        return;
    }

    n_qelscales = atoi(file.buffer[0][0]) / n_lambdas;
    q_elscales = (q_elscale_t*) malloc (n_qelscales * n_lambdas * sizeof(q_elscale_t));

    for (int i = 0; i < n_qelscales; i++) {
        for (int j = 0; j < n_lambdas; j++) {
            char *eptr;
            q_elscales[i + j * n_qelscales].qi = atoi(file.buffer[i + j * n_qelscales + 1][0]);
            q_elscales[i + j * n_qelscales].qj = atoi(file.buffer[i + j * n_qelscales + 1][1]);
            q_elscales[i + j * n_qelscales].mu = strtod(file.buffer[i + j * n_qelscales + 1][2], &eptr);
        }
    }

    clean_csv(file);
}

void init_qexclpairs(const char*filename) {
    csvfile_t file = read_csv(filename, 0, base_folder);

    if (file.n_lines < 1) {
        clean_csv(file);
        return;
    }

    n_qexclpairs = atoi(file.buffer[0][0]) / n_lambdas;
    q_exclpairs = (q_exclpair_t*) malloc (n_qexclpairs * n_lambdas * sizeof(q_exclpair_t));

    for (int i = 0; i < n_qexclpairs; i++) {
        for (int j = 0; j < n_lambdas; j++) {
            q_exclpairs[i + j * n_qexclpairs].ai = atoi(file.buffer[i + j * n_qexclpairs + 1][0]);
            q_exclpairs[i + j * n_qexclpairs].aj = atoi(file.buffer[i + j * n_qexclpairs + 1][1]);
            q_exclpairs[i + j * n_qexclpairs].excl = atoi(file.buffer[i + j * n_qexclpairs + 1][2]);
        }
    }

    clean_csv(file);
}

void init_qimpropers(const char*filename) {
    csvfile_t file = read_csv(filename, 0, base_folder);

    if (file.n_lines < 1) {
        clean_csv(file);
        return;
    }

    n_qimpropers = atoi(file.buffer[0][0]) / n_lambdas;
    q_impropers = (q_improper_t*) malloc (n_qimpropers * n_lambdas * sizeof(q_improper_t));

    for (int i = 0; i < n_qimpropers; i++) {
        for (int j = 0; j < n_lambdas; j++) {
            q_impropers[i + j * n_qimpropers].ai = atoi(file.buffer[i + j * n_qimpropers + 1][0]);
            q_impropers[i + j * n_qimpropers].aj = atoi(file.buffer[i + j * n_qimpropers + 1][1]);
            q_impropers[i + j * n_qimpropers].ak = atoi(file.buffer[i + j * n_qimpropers + 1][2]);
            q_impropers[i + j * n_qimpropers].al = atoi(file.buffer[i + j * n_qimpropers + 1][3]);
            q_impropers[i + j * n_qimpropers].code = atoi(file.buffer[i + j * n_qimpropers + 1][4]);
        }
    }

    clean_csv(file);
}

void init_qshakes(const char*filename) {
    csvfile_t file = read_csv(filename, 0, base_folder);

    if (file.n_lines < 1) {
        clean_csv(file);
        return;
    }

    n_qshakes = atoi(file.buffer[0][0]) / n_lambdas;
    q_shakes = (q_shake_t*) malloc (n_qshakes * n_lambdas * sizeof(q_shake_t));

    for (int i = 0; i < n_qshakes; i++) {
        for (int j = 0; j < n_lambdas; j++) {
            char *eptr;
            q_shakes[i + j * n_qshakes].ai = atoi(file.buffer[i + j * n_qshakes + 1][0]);
            q_shakes[i + j * n_qshakes].aj = atoi(file.buffer[i + j * n_qshakes + 1][1]);
            q_shakes[i + j * n_qshakes].dist = strtod(file.buffer[i + j * n_qshakes + 1][2], &eptr);
        }
    }

    clean_csv(file);
}

void init_qsoftcores(const char*filename) {
    csvfile_t file = read_csv(filename, 0, base_folder);

    if (file.n_lines < 1) {
        clean_csv(file);
        return;
    }

    printf("file.n_lines = %d\n", file.n_lines);

    n_qsoftcores = atoi(file.buffer[0][0]) / n_lambdas;
    q_softcores = (q_softcore_t*) malloc (n_qsoftcores * n_lambdas * sizeof(q_softcore_t));

    for (int i = 0; i < n_qsoftcores; i++) {
        for (int j = 0; j < n_lambdas; j++) {
            char *eptr;
            q_softcores[i + j * n_qsoftcores].s = strtod(file.buffer[i + j * n_qsoftcores + 1][0], &eptr);
        }
    }

    clean_csv(file);
}

void init_qtorsions(const char*filename) {
    csvfile_t file = read_csv(filename, 0, base_folder);

    if (file.n_lines < 1) {
        clean_csv(file);
        return;
    }

    n_qtorsions = atoi(file.buffer[0][0]) / n_lambdas;
    q_torsions = (q_torsion_t*) malloc (n_qtorsions * n_lambdas * sizeof(q_torsion_t));

    for (int i = 0; i < n_qtorsions; i++) {
        for (int j = 0; j < n_lambdas; j++) {
            q_torsions[i + j * n_qtorsions].ai = atoi(file.buffer[i + j * n_qtorsions + 1][0]);
            q_torsions[i + j * n_qtorsions].aj = atoi(file.buffer[i + j * n_qtorsions + 1][1]);
            q_torsions[i + j * n_qtorsions].ak = atoi(file.buffer[i + j * n_qtorsions + 1][2]);
            q_torsions[i + j * n_qtorsions].al = atoi(file.buffer[i + j * n_qtorsions + 1][3]);
            q_torsions[i + j * n_qtorsions].code = atoi(file.buffer[i + j * n_qtorsions + 1][4]);
        }
    }

    clean_csv(file);
}

/* =============================================
 * == FROM INPUT FILE
 * =============================================
 */

void init_icoords(const char*filename) {
    csvfile_t file = read_csv(filename, 0, base_folder);

    if (file.n_lines < 1) {
        clean_csv(file);
        return;
    }

    for (int i = 0; i < n_atoms; i++) {
        char *eptr;
        coords[i].x = strtod(file.buffer[i+1][0], &eptr);
        coords[i].y = strtod(file.buffer[i+1][1], &eptr);
        coords[i].z = strtod(file.buffer[i+1][2], &eptr);
    }

    clean_csv(file);
}

void init_ivelocities(const char*filename) {
    csvfile_t file = read_csv(filename, 0, base_folder);

    if (file.n_lines < 1) {
        clean_csv(file);
        return;
    }

    for (int i = 0; i < n_atoms; i++) {
        char *eptr;
        velocities[i].x = strtod(file.buffer[i+1][0], &eptr);
        velocities[i].y = strtod(file.buffer[i+1][1], &eptr);
        velocities[i].z = strtod(file.buffer[i+1][2], &eptr);
    }

    clean_csv(file);
}