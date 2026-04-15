#include "common/include/parse.h"

#include <stdio.h>
#include <unistd.h>

#include <cmath>
#include <cstring>

#include "common/include/constants.h"
#include "common/include/context.h"

csvfile_t read_csv(const char* filename, int ext, const char* base_folder) {
    csvfile_t retval;

    retval.ext = ext;

    char path[1024];
    sprintf(path, "%s/%s", base_folder, filename);
    if (access(path, F_OK) == -1) {
        printf(">>> FATAL: The following file could not be found. Exiting...\n");
        puts(path);
        exit(EXIT_FAILURE);
    }

    // File handle
    FILE* fp;

    fp = fopen(path, "r");
    if (fp == NULL)
        exit(EXIT_FAILURE);

    // Get number of lines
    char nlines[COLUMN_WIDTH];

    if (fgets(nlines, COLUMN_WIDTH, fp)) {
        retval.n_lines = atoi(nlines);
    } else {
        retval.n_lines = 0;
        return retval;
    }

    if (retval.n_lines == 0) {
        return retval;
    }

    char line[N_COLUMNS * COLUMN_WIDTH];
    retval.buffer = (char***)malloc(retval.n_lines * N_COLUMNS * sizeof(char**));

    for (int i = 0; i <= retval.n_lines + ext; i++) {
        retval.buffer[i] = (char**)malloc(N_COLUMNS * sizeof(char*));
        for (int j = 0; j < N_COLUMNS; j++) {
            retval.buffer[i][j] = (char*)malloc(COLUMN_WIDTH * sizeof(char));
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
             tok = strtok(NULL, ";\n")) {
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

void init_md(const char* filename) {
    auto& ctx = Context::instance();
    csvfile_t file = read_csv(filename, 0, ctx.base_folder.c_str());
    auto& md = ctx.md;
    char* eptr;

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
    ctx.n_lambdas = atoi(file.buffer[k][0]);
#ifdef VERBOSE
    printf("reading in %d lambdas (%s in file)\n", ctx.n_lambdas, file.buffer[k][1]);
#endif
    ctx.lambdas.resize(ctx.n_lambdas);
    k++;
    for (int i = 0; i < ctx.n_lambdas; i++) {
        ctx.lambdas[i] = strtod(file.buffer[k][0], &eptr);
        k++;
    }

    // [sequence_restraints]
    printf("k = %d\n", k);
    ctx.n_restrseqs = atoi(file.buffer[k][0]);
    printf("reading in %d sequence restraints (%s in file)\n", ctx.n_restrseqs, file.buffer[k][1]);
    ctx.restrseqs.resize(ctx.n_restrseqs);
    k++;
    for (int i = 0; i < ctx.n_restrseqs; i++) {
        restrseq_t restrseq;

        restrseq.ai = atoi(file.buffer[k][0]);
        restrseq.aj = atoi(file.buffer[k][1]);
        restrseq.k = strtod(file.buffer[k][2], &eptr);
        restrseq.ih = strcmp(file.buffer[k][3], "1") == 0;
        restrseq.to_center = atoi(file.buffer[k][4]);

        ctx.restrseqs[i] = restrseq;
        k++;
    }

    // [position_restraints]
    ctx.n_restrspos = atoi(file.buffer[k][0]);
    printf("reading in %d position restraints\n (%s in file )", ctx.n_restrspos, file.buffer[k][1]);
    ctx.restrspos.resize(ctx.n_restrspos);
    k++;
    for (int i = 0; i < ctx.n_restrspos; i++) {
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

        ctx.restrspos[i] = restrpos;
        k++;
    }

    // [distance_restraints]
    ctx.n_restrdists = atoi(file.buffer[k][0]);
    ctx.restrdists.resize(ctx.n_restrdists);
    printf("reading in %d distance restraints (%s in file)\n", ctx.n_restrdists, file.buffer[k][1]);
    k++;
    for (int i = 0; i < ctx.n_restrdists; i++) {
        restrdis_t restrdist;

        restrdist.ai = atoi(file.buffer[k][0]);
        restrdist.aj = atoi(file.buffer[k][1]);
        restrdist.d1 = strtod(file.buffer[k][2], &eptr);
        restrdist.d2 = strtod(file.buffer[k][3], &eptr);
        restrdist.k = strtod(file.buffer[k][4], &eptr);
        restrdist.ipsi = atoi(file.buffer[k][5]);

        ctx.restrdists[i] = restrdist;
        k++;
    }

    // [angle_restraints]
    ctx.n_restrangs = atoi(file.buffer[k][0]);
    ctx.restrangs.resize(ctx.n_restrangs);
    printf("reading in %d angle restraints (%s in file)\n", ctx.n_restrangs, file.buffer[k][1]);
    k++;
    for (int i = 0; i < ctx.n_restrangs; i++) {
        restrang_t restrang;

        restrang.ai = atoi(file.buffer[k][0]);
        restrang.aj = atoi(file.buffer[k][1]);
        restrang.ak = atoi(file.buffer[k][2]);
        restrang.ipsi = atoi(file.buffer[k][3]);
        restrang.ang = strtod(file.buffer[k][4], &eptr);
        restrang.k = strtod(file.buffer[k][5], &eptr);

        ctx.restrangs[i] = restrang;
        k++;
    }

    // [wall_restraints]
    ctx.n_restrwalls = atoi(file.buffer[k][0]);
    ctx.restrwalls.resize(ctx.n_restrwalls);
    printf("reading in %d wall restraints (%s in file)\n", ctx.n_restrwalls, file.buffer[k][1]);
    k++;
    for (int i = 0; i < ctx.n_restrwalls; i++) {
        restrwall_t restrwall;

        restrwall.ai = atoi(file.buffer[k][0]);
        restrwall.aj = atoi(file.buffer[k][1]);
        restrwall.d = strtod(file.buffer[k][2], &eptr);
        restrwall.k = strtod(file.buffer[k][3], &eptr);
        restrwall.dMorse = strtod(file.buffer[k][4], &eptr);
        restrwall.aMorse = strtod(file.buffer[k][5], &eptr);
        restrwall.ih = strcmp(file.buffer[k][6], "1") == 0;

        ctx.restrwalls[i] = restrwall;
        k++;
    }

    clean_csv(file);
}

/* =============================================
 * == FROM TOPOLOGY FILE
 * =============================================
 */

void init_topo(const char* filename) {
    auto& ctx = Context::instance();
    csvfile_t file = read_csv(filename, 0, ctx.base_folder.c_str());
    auto& topo = ctx.topo;
    char* eptr;

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
    topo.vdw_rule = atoi(file.buffer[8][0]);
    if (topo.vdw_rule == 0) {
        topo.vdw_rule = VDW_GEOMETRIC;
    }
#ifdef VERBOSE
    printf("read %d into vdw_rule (1=geometric, 2=arithmetic)\n", topo.vdw_rule);
#endif

    if (topo.vdw_rule != VDW_GEOMETRIC && topo.vdw_rule != VDW_ARITHMETIC) {
        printf(">>> FATAL: Invalid vdw_rule %d. Must be 1 (geometric) or 2 (arithmetic). Exiting...\n", topo.vdw_rule);
        exit(EXIT_FAILURE);
    }

    clean_csv(file);
}

void init_coords(const char* filename) {
    auto& ctx = Context::instance();
    csvfile_t file = read_csv(filename, 1, ctx.base_folder.c_str());

    ctx.n_coords = 0;
    ctx.n_atoms = 0;
    ctx.n_atoms_solute = 0;

    if (file.n_lines < 1) {
        clean_csv(file);
        return;
    }

    ctx.n_coords = atoi(file.buffer[0][0]);
    ctx.n_atoms = ctx.n_coords;

    ctx.n_atoms_solute = atoi(file.buffer[1][0]);

    ctx.coords.resize(ctx.n_atoms);
    ctx.coords_top.resize(ctx.n_atoms);

    for (int i = 0; i < file.n_lines; i++) {
        char* eptr;

        ctx.coords[i].x = strtod(file.buffer[i + 2][0], &eptr);
        ctx.coords[i].y = strtod(file.buffer[i + 2][1], &eptr);
        ctx.coords[i].z = strtod(file.buffer[i + 2][2], &eptr);

        ctx.coords_top[i] = ctx.coords[i];
    }

    clean_csv(file);
}

void init_bonds(const char* filename) {
    auto& ctx = Context::instance();
    csvfile_t file = read_csv(filename, 1, ctx.base_folder.c_str());

    ctx.n_bonds = 0;
    ctx.n_bonds_solute = 0;

    if (file.n_lines < 1) {
        clean_csv(file);
        return;
    }

    ctx.n_bonds = atoi(file.buffer[0][0]);
    ctx.n_bonds_solute = atoi(file.buffer[1][0]);

    ctx.bonds.resize(ctx.n_bonds);

    for (int i = 0; i < ctx.n_bonds; i++) {
        bond_t bond;

        bond.ai = atoi(file.buffer[i + 2][0]);
        bond.aj = atoi(file.buffer[i + 2][1]);
        bond.code = atoi(file.buffer[i + 2][2]);

        ctx.bonds[i] = bond;
    }

    clean_csv(file);
}

void init_cbonds(const char* filename) {
    auto& ctx = Context::instance();
    csvfile_t file = read_csv(filename, 0, ctx.base_folder.c_str());

    ctx.n_cbonds = 0;

    if (file.n_lines < 1) {
        clean_csv(file);
        return;
    }

    ctx.n_cbonds = atoi(file.buffer[0][0]);
    ctx.cbonds.resize(ctx.n_cbonds);

    for (int i = 0; i < ctx.n_cbonds; i++) {
        cbond_t cbond;
        char* eptr;

        cbond.code = atoi(file.buffer[i + 1][0]);
        cbond.kb = strtod(file.buffer[i + 1][1], &eptr);
        cbond.b0 = strtod(file.buffer[i + 1][2], &eptr);

        ctx.cbonds[i] = cbond;
    }

    clean_csv(file);
}

void init_angles(const char* filename) {
    auto& ctx = Context::instance();
    csvfile_t file = read_csv(filename, 1, ctx.base_folder.c_str());

    ctx.n_angles = 0;
    ctx.n_angles_solute = 0;

    if (file.n_lines < 1) {
        clean_csv(file);
        return;
    }

    ctx.n_angles = atoi(file.buffer[0][0]);
    ctx.n_angles_solute = atoi(file.buffer[1][0]);

    ctx.angles.resize(ctx.n_angles);

    for (int i = 0; i < ctx.n_angles; i++) {
        angle_t angle;

        angle.ai = atoi(file.buffer[i + 2][0]);
        angle.aj = atoi(file.buffer[i + 2][1]);
        angle.ak = atoi(file.buffer[i + 2][2]);
        angle.code = atoi(file.buffer[i + 2][3]);

        ctx.angles[i] = angle;
    }

    clean_csv(file);
}

void init_cangles(const char* filename) {
    auto& ctx = Context::instance();
    csvfile_t file = read_csv(filename, 0, ctx.base_folder.c_str());

    ctx.n_cangles = 0;

    if (file.n_lines < 1) {
        clean_csv(file);
        return;
    }

    ctx.n_cangles = atoi(file.buffer[0][0]);
    ctx.cangles.resize(ctx.n_cangles);

    for (int i = 0; i < ctx.n_cangles; i++) {
        cangle_t cangle;
        char* eptr;

        cangle.code = atoi(file.buffer[i + 1][0]);
        cangle.kth = strtod(file.buffer[i + 1][1], &eptr);
        cangle.th0 = strtod(file.buffer[i + 1][2], &eptr);

        ctx.cangles[i] = cangle;
    }

    clean_csv(file);
}

void init_excluded(const char* filename) {
    auto& ctx = Context::instance();
    ctx.excluded = std::make_unique<bool[]>(ctx.n_atoms);
    for (int i = 0; i < ctx.n_atoms; ++i) {
        ctx.excluded[i] = false;
    }
    ctx.n_excluded = 0;

    FILE* fp;

    char path[1024];
    sprintf(path, "%s/%s", ctx.base_folder.c_str(), filename);

    if (access(path, F_OK) == -1) {
        printf(">>> FATAL: The following file could not be found. Exiting...\n");
        puts(path);
        exit(EXIT_FAILURE);
    }

    fp = fopen(path, "r");
    if (fp == NULL)
        exit(EXIT_FAILURE);

    // --- dynamically read a full line ---
    char* line = NULL;
    size_t len = 0;
    ssize_t read_len = getline(&line, &len, fp);
    if (read_len == -1) {
        fprintf(stderr, "Error: failed to read line from %s\n", path);
        free(line);
        fclose(fp);
        exit(EXIT_FAILURE);
    }

    // --- parse ---
    ctx.n_excluded = 0;
    for (int i = 0; i < ctx.n_atoms && i < read_len; i++) {
        bool excl = (line[i] == '1');
        ctx.excluded[i] = excl;
        if (excl) ctx.n_excluded++;
    }

    printf("Number of excluded atoms: %d\n", ctx.n_excluded);

    free(line);
    fclose(fp);
}

void init_torsions(const char* filename) {
    auto& ctx = Context::instance();
    csvfile_t file = read_csv(filename, 1, ctx.base_folder.c_str());

    ctx.n_torsions = 0;
    ctx.n_torsions_solute = 0;

    if (file.n_lines < 1) {
        clean_csv(file);
        return;
    }

    ctx.n_torsions = atoi(file.buffer[0][0]);
    ctx.n_torsions_solute = atoi(file.buffer[1][0]);

    ctx.torsions.resize(ctx.n_torsions);

    for (int i = 0; i < ctx.n_torsions; i++) {
        torsion_t torsion;

        torsion.ai = atoi(file.buffer[i + 2][0]);
        torsion.aj = atoi(file.buffer[i + 2][1]);
        torsion.ak = atoi(file.buffer[i + 2][2]);
        torsion.al = atoi(file.buffer[i + 2][3]);
        torsion.code = atoi(file.buffer[i + 2][4]);

        ctx.torsions[i] = torsion;
    }

    clean_csv(file);
}

void init_ctorsions(const char* filename) {
    auto& ctx = Context::instance();
    csvfile_t file = read_csv(filename, 0, ctx.base_folder.c_str());

    ctx.n_ctorsions = 0;

    if (file.n_lines < 1) {
        clean_csv(file);
        return;
    }

    ctx.n_ctorsions = atoi(file.buffer[0][0]);

    ctx.ctorsions.resize(ctx.n_ctorsions);

    for (int i = 0; i < ctx.n_ctorsions; i++) {
        ctorsion_t ctorsion;
        char* eptr;

        ctorsion.code = atoi(file.buffer[i + 1][0]);
        ctorsion.k = strtod(file.buffer[i + 1][1], &eptr);
        ctorsion.n = strtod(file.buffer[i + 1][2], &eptr);
        ctorsion.d = strtod(file.buffer[i + 1][3], &eptr);
        ctorsion.paths = strtod(file.buffer[i + 1][4], &eptr);
        ctorsion.paths = 1.0 / (ctorsion.paths);

        ctx.ctorsions[i] = ctorsion;
    }

    clean_csv(file);
}

void init_impropers(const char* filename) {
    auto& ctx = Context::instance();
    csvfile_t file = read_csv(filename, 1, ctx.base_folder.c_str());

    ctx.n_impropers = 0;
    ctx.n_impropers_solute = 0;

    if (file.n_lines < 1) {
        clean_csv(file);
        return;
    }

    ctx.n_impropers = atoi(file.buffer[0][0]);
    ctx.n_impropers_solute = atoi(file.buffer[1][0]);

    ctx.impropers.resize(ctx.n_impropers);

    for (int i = 0; i < ctx.n_impropers; i++) {
        improper_t improper;

        improper.ai = atoi(file.buffer[i + 2][0]);
        improper.aj = atoi(file.buffer[i + 2][1]);
        improper.ak = atoi(file.buffer[i + 2][2]);
        improper.al = atoi(file.buffer[i + 2][3]);
        improper.code = atoi(file.buffer[i + 2][4]);

        ctx.impropers[i] = improper;
    }

    clean_csv(file);
}

void init_cimpropers(const char* filename) {
    auto& ctx = Context::instance();
    csvfile_t file = read_csv(filename, 0, ctx.base_folder.c_str());

    ctx.n_cimpropers = 0;

    if (file.n_lines < 1) {
        clean_csv(file);
        return;
    }

    ctx.n_cimpropers = atoi(file.buffer[0][0]);

    ctx.cimpropers.resize(ctx.n_cimpropers);

    for (int i = 0; i < ctx.n_cimpropers; i++) {
        cimproper_t cimproper;
        char* eptr;

        cimproper.code = atoi(file.buffer[i + 1][0]);
        cimproper.k = strtod(file.buffer[i + 1][1], &eptr);
        cimproper.phi0 = strtod(file.buffer[i + 1][2], &eptr);

        ctx.cimpropers[i] = cimproper;
    }

    clean_csv(file);
}

void init_charges(const char* filename) {
    auto& ctx = Context::instance();
    csvfile_t file = read_csv(filename, 0, ctx.base_folder.c_str());

    ctx.n_charges = 0;

    if (file.n_lines < 1) {
        clean_csv(file);
        return;
    }

    ctx.n_charges = atoi(file.buffer[0][0]);

    ctx.charges.resize(ctx.n_charges);

    for (int i = 0; i < ctx.n_charges; i++) {
        charge_t charge;

        charge.a = atoi(file.buffer[i + 1][0]);
        charge.code = atoi(file.buffer[i + 1][1]);

        ctx.charges[i] = charge;
    }

    clean_csv(file);
}

void init_ccharges(const char* filename) {
    auto& ctx = Context::instance();
    csvfile_t file = read_csv(filename, 0, ctx.base_folder.c_str());

    ctx.n_ccharges = 0;

    if (file.n_lines < 1) {
        clean_csv(file);
        return;
    }

    ctx.n_ccharges = atoi(file.buffer[0][0]);

    ctx.ccharges.resize(ctx.n_ccharges);

    for (int i = 0; i < ctx.n_ccharges; i++) {
        ccharge_t ccharge;
        char* eptr;

        ccharge.code = atoi(file.buffer[i + 1][0]);
        ccharge.charge = strtod(file.buffer[i + 1][1], &eptr);

        ctx.ccharges[i] = ccharge;
    }

    clean_csv(file);
}

void init_ngbrs14(const char* filename) {
    auto& ctx = Context::instance();
    FILE* fp;

    char path[1024];
    sprintf(path, "%s/%s", ctx.base_folder.c_str(), filename);

    if (access(path, F_OK) == -1) {
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
    } else {
        return;
    }

    int lineI = 0;

    while (fgets(line, 1024, fp)) {
        for (int i = 0; i < line_width; i++) {
            if (line[i] == '1') {
                int ix = lineI;
                int jx = (lineI + i + 1) % lines;
                // if (ix < 100 && jx < 100) printf("i = %d j = %d\n", ix+1, jx+1);
                ctx.LJ_matrix[ix * ctx.n_atoms_solute + jx] = 1;
                ctx.LJ_matrix[jx * ctx.n_atoms_solute + ix] = 1;

                ctx.ngbrs_14.push_back({ix, jx, NONBONDED_14_PP}); // the type may is wrong, just set in here
            }
        }
        lineI++;
    }

    fclose(fp);
}

void init_ngbrs23(const char* filename) {
    auto& ctx = Context::instance();
    FILE* fp;

    char path[1024];
    sprintf(path, "%s/%s", ctx.base_folder.c_str(), filename);

    if (access(path, F_OK) == -1) {
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
    } else {
        return;
    }

    int lineI = 0;

    while (fgets(line, 1024, fp)) {
        for (int i = 0; i < line_width; i++) {
            if (line[i] == '1') {
                int ix = lineI;
                int jx = (lineI + i + 1) % lines;
                // if (ix < 100 && jx < 100) printf("i = %d j = %d\n", ix+1, jx+1);
                ctx.LJ_matrix[ix * ctx.n_atoms_solute + jx] = 3;
                ctx.LJ_matrix[jx * ctx.n_atoms_solute + ix] = 3;
            }
        }
        lineI++;
    }

    fclose(fp);
}

void init_ngbrs14_long(const char* filename) {
    auto& ctx = Context::instance();
    csvfile_t file = read_csv(filename, 0, ctx.base_folder.c_str());

    if (file.n_lines < 1) {
        clean_csv(file);
        return;
    }

    int n_ngbrs14_long = atoi(file.buffer[0][0]);

    for (int i = 0; i < n_ngbrs14_long; i++) {
        int ix = atoi(file.buffer[i + 1][0]) - 1;
        int jx = atoi(file.buffer[i + 1][1]) - 1;
        ctx.LJ_matrix[ix * ctx.n_atoms_solute + jx] = 1;
        ctx.LJ_matrix[jx * ctx.n_atoms_solute + ix] = 1;


        int mi_x = std::min(ix, jx);
        int mx_x = std::max(ix, jx);
        ctx.ngbrs_14.push_back({mi_x, mx_x, NONBONDED_14_PP}); // the type may is wrong, just set in here

    }

    clean_csv(file);
}

void init_ngbrs23_long(const char* filename) {
    auto& ctx = Context::instance();
    csvfile_t file = read_csv(filename, 0, ctx.base_folder.c_str());

    if (file.n_lines < 1) {
        clean_csv(file);
        return;
    }

    int n_ngbrs23_long = atoi(file.buffer[0][0]);

    for (int i = 0; i < n_ngbrs23_long; i++) {
        int ix = atoi(file.buffer[i + 1][0]) - 1;
        int jx = atoi(file.buffer[i + 1][1]) - 1;
        ctx.LJ_matrix[ix * ctx.n_atoms_solute + jx] = 3;
        ctx.LJ_matrix[jx * ctx.n_atoms_solute + ix] = 3;
    }

    clean_csv(file);
}

void init_LJ_matrix() {
    auto& ctx = Context::instance();
    ctx.LJ_matrix.assign(ctx.n_atoms_solute * ctx.n_atoms_solute, 0);
}

void init_catypes(const char* filename) {
    auto& ctx = Context::instance();
    csvfile_t file = read_csv(filename, 0, ctx.base_folder.c_str());

    ctx.n_catypes = 0;

    if (file.n_lines < 1) {
        clean_csv(file);
        return;
    }

    ctx.n_catypes = atoi(file.buffer[0][0]);
    ctx.catypes.resize(ctx.n_catypes);

    for (int i = 0; i < ctx.n_catypes; i++) {
        catype_t catype;
        char* eptr;

        catype.code = atoi(file.buffer[i + 1][0]);
        catype.m = strtod(file.buffer[i + 1][1], &eptr);
        catype.aii_normal = strtod(file.buffer[i + 1][2], &eptr);
        catype.bii_normal = strtod(file.buffer[i + 1][3], &eptr);
        // Legacy polar columns are currently unused.
        strtod(file.buffer[i + 1][4], &eptr);
        strtod(file.buffer[i + 1][5], &eptr);
        catype.aii_1_4 = strtod(file.buffer[i + 1][6], &eptr);
        catype.bii_1_4 = strtod(file.buffer[i + 1][7], &eptr);

        ctx.catypes[i] = catype;
    }

    // Preprocess bii parameters for arithmetic rule: convert ε to √ε
    // This matches Fortran md.f90:14747-14752 preprocessing
    if (ctx.topo.vdw_rule == VDW_ARITHMETIC) {
        for (int i = 0; i < ctx.n_catypes; i++) {
            ctx.catypes[i].bii_normal = sqrt(fabs(ctx.catypes[i].bii_normal));
            ctx.catypes[i].bii_1_4 = sqrt(fabs(ctx.catypes[i].bii_1_4));
        }
#ifdef VERBOSE
        printf("Preprocessed catypes bii parameters for arithmetic vdW rule\n");
#endif
    }

    clean_csv(file);
}

void init_atypes(const char* filename) {
    auto& ctx = Context::instance();
    csvfile_t file = read_csv(filename, 0, ctx.base_folder.c_str());

    ctx.n_atypes = 0;

    if (file.n_lines < 1) {
        clean_csv(file);
        return;
    }

    ctx.n_atypes = atoi(file.buffer[0][0]);

    ctx.atypes.resize(ctx.n_atypes);
    for (int i = 0; i < ctx.n_atypes; i++) {
        atype_t atype;

        atype.a = atoi(file.buffer[i + 1][0]);
        atype.code = atoi(file.buffer[i + 1][1]);

        ctx.atypes[i] = atype;
    }

    clean_csv(file);
}

void init_molecules(const char* filename) {
    auto& ctx = Context::instance();
    csvfile_t file = read_csv(filename, 0, ctx.base_folder.c_str());

    ctx.n_molecules = 0;

    if (file.n_lines < 1) {
        clean_csv(file);
        return;
    }

    ctx.n_molecules = atoi(file.buffer[0][0]);
    ctx.molecules.resize(ctx.n_molecules);

    for (int i = 0; i < ctx.n_molecules; i++) {
        ctx.molecules[i] = atoi(file.buffer[i + 1][0]);
    }

    clean_csv(file);
}

void init_charge_groups(const char* filename) {
    auto& ctx = Context::instance();
    csvfile_t file = read_csv(filename, 0, ctx.base_folder.c_str());

    if (file.n_lines < 1) {
        clean_csv(file);
        return;
    }

    ctx.n_cgrps_solute = atoi(file.buffer[1][0]);
    ctx.n_cgrps_solvent = atoi(file.buffer[1][1]);
    ctx.iuse_switch_atom = atoi(file.buffer[1][2]);

    int n_charge_groups = ctx.n_cgrps_solute + ctx.n_cgrps_solvent;

    ctx.charge_groups.resize(n_charge_groups);

    int line_nr = 2;
    int n_atoms_crgp = 0;

    for (int i = 0; i < n_charge_groups; i++) {
        cgrp_t charge_group;

        n_atoms_crgp = atoi(file.buffer[line_nr][0]);
        charge_group.n_atoms = n_atoms_crgp;
        charge_group.iswitch = atoi(file.buffer[line_nr][1]);
        charge_group.a = new int[n_atoms_crgp];

        line_nr++;
        for (int j = 0; j < charge_group.n_atoms; j++) {
            charge_group.a[j] = atoi(file.buffer[line_nr][0]);
            line_nr++;
        }

        ctx.charge_groups[i] = charge_group;
    }

    clean_csv(file);
}

/* =============================================
 * == FROM FEP FILE
 * =============================================
 */

void init_qangcouples(const char* filename) {
    auto& ctx = Context::instance();
    csvfile_t file = read_csv(filename, 0, ctx.base_folder.c_str());

    ctx.n_qangcouples = 0;

    if (file.n_lines < 1) {
        clean_csv(file);
        return;
    }

    ctx.n_qangcouples = atoi(file.buffer[0][0]);
    ctx.q_angcouples.resize(ctx.n_qangcouples);

    for (int i = 0; i < ctx.n_qangcouples; i++) {
        ctx.q_angcouples[i].acode = atoi(file.buffer[i + 1][0]);
        ctx.q_angcouples[i].bcode = atoi(file.buffer[i + 1][1]);
    }

    clean_csv(file);
}

void init_qatoms(const char* filename) {
    auto& ctx = Context::instance();
    csvfile_t file = read_csv(filename, 0, ctx.base_folder.c_str());

    ctx.n_qatoms = 0;
    ctx.q_atoms.clear();

    if (file.n_lines < 1) {
        clean_csv(file);
        return;
    }

    ctx.n_qatoms = atoi(file.buffer[0][0]);
    ctx.q_atoms.resize(ctx.n_qatoms);

    for (int i = 0; i < ctx.n_qatoms; i++) {
        ctx.q_atoms[i] = atoi(file.buffer[i + 1][0]) - 1;
    }

    clean_csv(file);
}

void init_qcangles(const char* filename) {
    auto& ctx = Context::instance();
    csvfile_t file = read_csv(filename, 0, ctx.base_folder.c_str());

    ctx.n_qcangles = 0;

    if (file.n_lines < 1) {
        clean_csv(file);
        return;
    }

    ctx.n_qcangles = atoi(file.buffer[0][0]);
    ctx.q_cangles.resize(ctx.n_qcangles);

    for (int i = 0; i < ctx.n_qcangles; i++) {
        char* eptr;
        ctx.q_cangles[i].kth = strtod(file.buffer[i + 1][0], &eptr);
        ctx.q_cangles[i].th0 = strtod(file.buffer[i + 1][1], &eptr);
    }

    clean_csv(file);
}

void init_qcatypes(const char* filename) {
    auto& ctx = Context::instance();
    csvfile_t file = read_csv(filename, 0, ctx.base_folder.c_str());

    ctx.n_qcatypes = 0;

    if (file.n_lines < 1) {
        clean_csv(file);
        return;
    }

    ctx.n_qcatypes = atoi(file.buffer[0][0]);
    ctx.q_catypes.resize(ctx.n_qcatypes);

    for (int i = 0; i < ctx.n_qcatypes; i++) {
        char* eptr;
        ctx.q_catypes[i].code = i + 1;
        ctx.q_catypes[i].aii_normal = strtod(file.buffer[i + 1][1], &eptr);
        ctx.q_catypes[i].bii_normal = strtod(file.buffer[i + 1][2], &eptr);
        // Legacy q_catype Ci/ai columns are currently unused.
        strtod(file.buffer[i + 1][3], &eptr);
        strtod(file.buffer[i + 1][4], &eptr);
        ctx.q_catypes[i].aii_1_4 = strtod(file.buffer[i + 1][5], &eptr);
        ctx.q_catypes[i].bii_1_4 = strtod(file.buffer[i + 1][6], &eptr);
        ctx.q_catypes[i].m = strtod(file.buffer[i + 1][7], &eptr);
    }

    // Preprocess Bi parameters for arithmetic rule: convert ε to √ε
    if (ctx.topo.vdw_rule == VDW_ARITHMETIC) {
        for (int i = 0; i < ctx.n_qcatypes; i++) {
            ctx.q_catypes[i].bii_normal = sqrt(fabs(ctx.q_catypes[i].bii_normal));
            ctx.q_catypes[i].bii_1_4 = sqrt(fabs(ctx.q_catypes[i].bii_1_4));
        }
#ifdef VERBOSE
        printf("Preprocessed q_catypes Bi parameters for arithmetic vdW rule\n");
#endif
    }

    clean_csv(file);
}

void init_qcbonds(const char* filename) {
    auto& ctx = Context::instance();
    csvfile_t file = read_csv(filename, 0, ctx.base_folder.c_str());

    ctx.n_qcbonds = 0;

    if (file.n_lines < 1) {
        clean_csv(file);
        return;
    }

    ctx.n_qcbonds = atoi(file.buffer[0][0]);
    ctx.q_cbonds.resize(ctx.n_qcbonds);

    for (int i = 0; i < ctx.n_qcbonds; i++) {
        char* eptr;
        ctx.q_cbonds[i].kb = strtod(file.buffer[i + 1][0], &eptr);
        ctx.q_cbonds[i].b0 = strtod(file.buffer[i + 1][1], &eptr);
    }

    clean_csv(file);
}

void init_qcimpropers(const char* filename) {
    auto& ctx = Context::instance();
    csvfile_t file = read_csv(filename, 0, ctx.base_folder.c_str());

    ctx.n_qcimpropers = 0;

    if (file.n_lines < 1) {
        clean_csv(file);
        return;
    }

    ctx.n_qcimpropers = atoi(file.buffer[0][0]);
    ctx.q_cimpropers.resize(ctx.n_qcimpropers);

    for (int i = 0; i < ctx.n_qcimpropers; i++) {
        char* eptr;
        ctx.q_cimpropers[i].k = strtod(file.buffer[i + 1][0], &eptr);
        ctx.q_cimpropers[i].phi0 = strtod(file.buffer[i + 1][1], &eptr);
    }

    clean_csv(file);
}

void init_qctorsions(const char* filename) {
    auto& ctx = Context::instance();
    csvfile_t file = read_csv(filename, 0, ctx.base_folder.c_str());

    ctx.n_qctorsions = 0;

    if (file.n_lines < 1) {
        clean_csv(file);
        return;
    }

    ctx.n_qctorsions = atoi(file.buffer[0][0]);
    ctx.q_ctorsions.resize(ctx.n_qctorsions);

    for (int i = 0; i < ctx.n_qctorsions; i++) {
        char* eptr;
        ctx.q_ctorsions[i].k = strtod(file.buffer[i + 1][0], &eptr);
        ctx.q_ctorsions[i].n = strtod(file.buffer[i + 1][1], &eptr);
        ctx.q_ctorsions[i].d = strtod(file.buffer[i + 1][2], &eptr);
    }

    clean_csv(file);
}

void init_qoffdiags(const char* filename) {
    auto& ctx = Context::instance();
    csvfile_t file = read_csv(filename, 0, ctx.base_folder.c_str());

    ctx.n_qoffdiags = 0;

    if (file.n_lines < 1) {
        clean_csv(file);
        return;
    }

    ctx.n_qoffdiags = atoi(file.buffer[0][0]);
    ctx.q_offdiags.resize(ctx.n_qoffdiags);

    for (int i = 0; i < ctx.n_qoffdiags; i++) {
        char* eptr;
        ctx.q_offdiags[i].i = atoi(file.buffer[i + 1][0]);
        ctx.q_offdiags[i].j = atoi(file.buffer[i + 1][1]);
        ctx.q_offdiags[i].qk = atoi(file.buffer[i + 1][2]);
        ctx.q_offdiags[i].ql = atoi(file.buffer[i + 1][3]);
        ctx.q_offdiags[i].Aij = strtod(file.buffer[i + 1][4], &eptr);
        ctx.q_offdiags[i].muij = strtod(file.buffer[i + 1][5], &eptr);
    }

    clean_csv(file);
}

void init_qimprcouples(const char* filename) {
    auto& ctx = Context::instance();
    csvfile_t file = read_csv(filename, 0, ctx.base_folder.c_str());

    ctx.n_qimprcouples = 0;

    if (file.n_lines < 1) {
        clean_csv(file);
        return;
    }

    ctx.n_qimprcouples = atoi(file.buffer[0][0]);
    ctx.q_imprcouples.resize(ctx.n_qimprcouples);

    for (int i = 0; i < ctx.n_qimprcouples; i++) {
        ctx.q_imprcouples[i].icode = atoi(file.buffer[i + 1][0]);
        ctx.q_imprcouples[i].bcode = atoi(file.buffer[i + 1][1]);
    }

    clean_csv(file);
}

void init_qsoftpairs(const char* filename) {
    auto& ctx = Context::instance();
    csvfile_t file = read_csv(filename, 0, ctx.base_folder.c_str());

    ctx.n_qsoftpairs = 0;

    if (file.n_lines < 1) {
        clean_csv(file);
        return;
    }

    ctx.n_qsoftpairs = atoi(file.buffer[0][0]);
    ctx.q_softpairs.resize(ctx.n_qsoftpairs);

    for (int i = 0; i < ctx.n_qsoftpairs; i++) {
        ctx.q_softpairs[i].qi = atoi(file.buffer[i + 1][0]);
        ctx.q_softpairs[i].qj = atoi(file.buffer[i + 1][1]);
    }

    clean_csv(file);
}

void init_qtorcouples(const char* filename) {
    auto& ctx = Context::instance();
    csvfile_t file = read_csv(filename, 0, ctx.base_folder.c_str());

    ctx.n_qtorcouples = 0;

    if (file.n_lines < 1) {
        clean_csv(file);
        return;
    }

    ctx.n_qtorcouples = atoi(file.buffer[0][0]);
    ctx.q_torcouples.resize(ctx.n_qtorcouples);

    for (int i = 0; i < ctx.n_qtorcouples; i++) {
        ctx.q_torcouples[i].tcode = atoi(file.buffer[i + 1][0]);
        ctx.q_torcouples[i].bcode = atoi(file.buffer[i + 1][1]);
    }

    clean_csv(file);
}

void init_qangles(const char* filename) {
    auto& ctx = Context::instance();
    csvfile_t file = read_csv(filename, 0, ctx.base_folder.c_str());

    if (file.n_lines < 1) {
        clean_csv(file);
        return;
    }

    ctx.n_qangles = atoi(file.buffer[0][0]) / ctx.n_lambdas;
    ctx.q_angles.resize(ctx.n_qangles * ctx.n_lambdas);

    for (int i = 0; i < ctx.n_qangles; i++) {
        for (int j = 0; j < ctx.n_lambdas; j++) {
            ctx.q_angles[i + j * ctx.n_qangles].ai = atoi(file.buffer[i + j * ctx.n_qangles + 1][0]);
            ctx.q_angles[i + j * ctx.n_qangles].aj = atoi(file.buffer[i + j * ctx.n_qangles + 1][1]);
            ctx.q_angles[i + j * ctx.n_qangles].ak = atoi(file.buffer[i + j * ctx.n_qangles + 1][2]);
            ctx.q_angles[i + j * ctx.n_qangles].code = atoi(file.buffer[i + j * ctx.n_qangles + 1][3]);
        }
    }

    clean_csv(file);
}

void init_qatypes(const char* filename) {
    auto& ctx = Context::instance();
    csvfile_t file = read_csv(filename, 0, ctx.base_folder.c_str());

    if (file.n_lines < 1) {
        clean_csv(file);
        return;
    }

    ctx.q_atypes.resize(ctx.n_qatoms * ctx.n_lambdas);
    for (int i = 0; i < ctx.n_qatoms; i++) {
        for (int j = 0; j < ctx.n_lambdas; j++) {
            ctx.q_atypes[i + j * ctx.n_qatoms].code = atoi(file.buffer[i + j * ctx.n_qatoms + 1][0]);
        }
    }

    clean_csv(file);
}

void init_qbonds(const char* filename) {
    auto& ctx = Context::instance();
    csvfile_t file = read_csv(filename, 0, ctx.base_folder.c_str());

    if (file.n_lines < 1) {
        clean_csv(file);
        return;
    }

    ctx.n_qbonds = atoi(file.buffer[0][0]) / ctx.n_lambdas;
    ctx.q_bonds.resize(ctx.n_qbonds * ctx.n_lambdas);

    for (int i = 0; i < ctx.n_qbonds; i++) {
        for (int j = 0; j < ctx.n_lambdas; j++) {
            ctx.q_bonds[i + j * ctx.n_qbonds].ai = atoi(file.buffer[i + j * ctx.n_qbonds + 1][0]);
            ctx.q_bonds[i + j * ctx.n_qbonds].aj = atoi(file.buffer[i + j * ctx.n_qbonds + 1][1]);
            ctx.q_bonds[i + j * ctx.n_qbonds].code = atoi(file.buffer[i + j * ctx.n_qbonds + 1][2]);
        }
    }

    clean_csv(file);
}

void init_qcharges(const char* filename) {
    auto& ctx = Context::instance();
    csvfile_t file = read_csv(filename, 0, ctx.base_folder.c_str());

    if (file.n_lines < 1) {
        clean_csv(file);
        return;
    }

    ctx.q_charges.resize(ctx.n_qatoms * ctx.n_lambdas);

    for (int i = 0; i < ctx.n_qatoms; i++) {
        for (int j = 0; j < ctx.n_lambdas; j++) {
            char* eptr;
            ctx.q_charges[i + j * ctx.n_qatoms].charge = strtod(file.buffer[i + j * ctx.n_qatoms + 1][0], &eptr);
        }
    }

    clean_csv(file);
}

void init_qelscales(const char* filename) {
    auto& ctx = Context::instance();
    csvfile_t file = read_csv(filename, 0, ctx.base_folder.c_str());

    if (file.n_lines < 1) {
        clean_csv(file);
        return;
    }

    ctx.n_qelscales = atoi(file.buffer[0][0]) / ctx.n_lambdas;
    ctx.q_elscales.resize(ctx.n_qelscales * ctx.n_lambdas);

    for (int i = 0; i < ctx.n_qelscales; i++) {
        for (int j = 0; j < ctx.n_lambdas; j++) {
            char* eptr;
            ctx.q_elscales[i + j * ctx.n_qelscales].qi = atoi(file.buffer[i + j * ctx.n_qelscales + 1][0]);
            ctx.q_elscales[i + j * ctx.n_qelscales].qj = atoi(file.buffer[i + j * ctx.n_qelscales + 1][1]);
            ctx.q_elscales[i + j * ctx.n_qelscales].mu = strtod(file.buffer[i + j * ctx.n_qelscales + 1][2], &eptr);
        }
    }

    clean_csv(file);
}

void init_qexclpairs(const char* filename) {
    auto& ctx = Context::instance();
    csvfile_t file = read_csv(filename, 0, ctx.base_folder.c_str());

    if (file.n_lines < 1) {
        clean_csv(file);
        return;
    }

    ctx.n_qexclpairs = atoi(file.buffer[0][0]) / ctx.n_lambdas;
    ctx.q_exclpairs.resize(ctx.n_qexclpairs * ctx.n_lambdas);

    for (int i = 0; i < ctx.n_qexclpairs; i++) {
        for (int j = 0; j < ctx.n_lambdas; j++) {
            ctx.q_exclpairs[i + j * ctx.n_qexclpairs].ai = atoi(file.buffer[i + j * ctx.n_qexclpairs + 1][0]);
            ctx.q_exclpairs[i + j * ctx.n_qexclpairs].aj = atoi(file.buffer[i + j * ctx.n_qexclpairs + 1][1]);
            ctx.q_exclpairs[i + j * ctx.n_qexclpairs].excl = atoi(file.buffer[i + j * ctx.n_qexclpairs + 1][2]);
        }
    }

    clean_csv(file);
}

void init_qimpropers(const char* filename) {
    auto& ctx = Context::instance();
    csvfile_t file = read_csv(filename, 0, ctx.base_folder.c_str());

    if (file.n_lines < 1) {
        clean_csv(file);
        return;
    }

    ctx.n_qimpropers = atoi(file.buffer[0][0]) / ctx.n_lambdas;
    ctx.q_impropers.resize(ctx.n_qimpropers * ctx.n_lambdas);

    for (int i = 0; i < ctx.n_qimpropers; i++) {
        for (int j = 0; j < ctx.n_lambdas; j++) {
            ctx.q_impropers[i + j * ctx.n_qimpropers].ai = atoi(file.buffer[i + j * ctx.n_qimpropers + 1][0]);
            ctx.q_impropers[i + j * ctx.n_qimpropers].aj = atoi(file.buffer[i + j * ctx.n_qimpropers + 1][1]);
            ctx.q_impropers[i + j * ctx.n_qimpropers].ak = atoi(file.buffer[i + j * ctx.n_qimpropers + 1][2]);
            ctx.q_impropers[i + j * ctx.n_qimpropers].al = atoi(file.buffer[i + j * ctx.n_qimpropers + 1][3]);
            ctx.q_impropers[i + j * ctx.n_qimpropers].code = atoi(file.buffer[i + j * ctx.n_qimpropers + 1][4]);
        }
    }

    clean_csv(file);
}

void init_qshakes(const char* filename) {
    auto& ctx = Context::instance();
    csvfile_t file = read_csv(filename, 0, ctx.base_folder.c_str());

    if (file.n_lines < 1) {
        clean_csv(file);
        return;
    }

    ctx.n_qshakes = atoi(file.buffer[0][0]) / ctx.n_lambdas;
    ctx.q_shakes.resize(ctx.n_qshakes * ctx.n_lambdas);

    for (int i = 0; i < ctx.n_qshakes; i++) {
        for (int j = 0; j < ctx.n_lambdas; j++) {
            char* eptr;
            ctx.q_shakes[i + j * ctx.n_qshakes].ai = atoi(file.buffer[i + j * ctx.n_qshakes + 1][0]);
            ctx.q_shakes[i + j * ctx.n_qshakes].aj = atoi(file.buffer[i + j * ctx.n_qshakes + 1][1]);
            ctx.q_shakes[i + j * ctx.n_qshakes].dist = strtod(file.buffer[i + j * ctx.n_qshakes + 1][2], &eptr);
        }
    }

    clean_csv(file);
}

void init_qsoftcores(const char* filename) {
    auto& ctx = Context::instance();
    csvfile_t file = read_csv(filename, 0, ctx.base_folder.c_str());

    if (file.n_lines < 1) {
        clean_csv(file);
        return;
    }

    printf("file.n_lines = %d\n", file.n_lines);

    ctx.n_qsoftcores = atoi(file.buffer[0][0]) / ctx.n_lambdas;
    ctx.q_softcores.resize(ctx.n_qsoftcores * ctx.n_lambdas);

    for (int i = 0; i < ctx.n_qsoftcores; i++) {
        for (int j = 0; j < ctx.n_lambdas; j++) {
            char* eptr;
            ctx.q_softcores[i + j * ctx.n_qsoftcores].s = strtod(file.buffer[i + j * ctx.n_qsoftcores + 1][0], &eptr);
        }
    }

    clean_csv(file);
}

void init_qtorsions(const char* filename) {
    auto& ctx = Context::instance();
    csvfile_t file = read_csv(filename, 0, ctx.base_folder.c_str());

    if (file.n_lines < 1) {
        clean_csv(file);
        return;
    }

    ctx.n_qtorsions = atoi(file.buffer[0][0]) / ctx.n_lambdas;
    ctx.q_torsions.resize(ctx.n_qtorsions * ctx.n_lambdas);

    for (int i = 0; i < ctx.n_qtorsions; i++) {
        for (int j = 0; j < ctx.n_lambdas; j++) {
            ctx.q_torsions[i + j * ctx.n_qtorsions].ai = atoi(file.buffer[i + j * ctx.n_qtorsions + 1][0]);
            ctx.q_torsions[i + j * ctx.n_qtorsions].aj = atoi(file.buffer[i + j * ctx.n_qtorsions + 1][1]);
            ctx.q_torsions[i + j * ctx.n_qtorsions].ak = atoi(file.buffer[i + j * ctx.n_qtorsions + 1][2]);
            ctx.q_torsions[i + j * ctx.n_qtorsions].al = atoi(file.buffer[i + j * ctx.n_qtorsions + 1][3]);
            ctx.q_torsions[i + j * ctx.n_qtorsions].code = atoi(file.buffer[i + j * ctx.n_qtorsions + 1][4]);
        }
    }

    clean_csv(file);
}

/* =============================================
 * == FROM INPUT FILE
 * =============================================
 */

void init_icoords(const char* filename) {
    auto& ctx = Context::instance();
    csvfile_t file = read_csv(filename, 0, ctx.base_folder.c_str());

    if (file.n_lines < 1) {
        clean_csv(file);
        return;
    }

    for (int i = 0; i < ctx.n_atoms; i++) {
        char* eptr;
        ctx.coords[i].x = strtod(file.buffer[i + 1][0], &eptr);
        ctx.coords[i].y = strtod(file.buffer[i + 1][1], &eptr);
        ctx.coords[i].z = strtod(file.buffer[i + 1][2], &eptr);
    }

    clean_csv(file);
}

void init_ivelocities(const char* filename) {
    auto& ctx = Context::instance();
    csvfile_t file = read_csv(filename, 0, ctx.base_folder.c_str());

    if (file.n_lines < 1) {
        clean_csv(file);
        return;
    }

    for (int i = 0; i < ctx.n_atoms; i++) {
        char* eptr;
        ctx.velocities[i].x = strtod(file.buffer[i + 1][0], &eptr);
        ctx.velocities[i].y = strtod(file.buffer[i + 1][1], &eptr);
        ctx.velocities[i].z = strtod(file.buffer[i + 1][2], &eptr);
    }

    clean_csv(file);
}
