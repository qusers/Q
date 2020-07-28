#include "parse.h"
#include "system.h"
#include <stdio.h>
#include <unistd.h>

csvfile_t read_csv(char* filename, int ext, char* base_folder) {
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

    char line[N_COLUMNS * COLUMN_WIDTH];
    retval.buffer = (char***) malloc(retval.n_lines * N_COLUMNS * sizeof(char**));

    for (int i = 0; i <= retval.n_lines + ext; i++) {
        retval.buffer[i] = (char**) malloc(N_COLUMNS * sizeof(char*));
        for (int j = 0; j < N_COLUMNS; j++) {
            retval.buffer[i][j] = (char*) malloc(N_COLUMNS * sizeof(char));
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
    for (int i = 0; i <= file.n_lines + file.ext; i++) {
        for (int j = 0; j < N_COLUMNS; j++) {
            free(file.buffer[i][j]);
        }
        free(file.buffer[i]);
    }
    free(file.buffer);
}


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