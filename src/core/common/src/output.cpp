#include "output.h"
#include "context.h"

void print_energies() {
    auto& host = Context::instance();
    auto *lambdas = host.lambdas->cpu_data_p;
    printf("[temperature]\n");
    printf("Temp\t%f\n", host.Temp);
    printf("\n");

    printf("[bonded]\n");
    printf("type\tUbond\tUangle\tUtor\tUimp\n");
    printf("p\t%f\t%f\t%f\t%f\n", host.E_bond_p.Ubond, host.E_bond_p.Uangle, host.E_bond_p.Utor, host.E_bond_p.Uimp);
    printf("w\t%f\t%f\t%f\t%f\n", host.E_bond_w.Ubond, host.E_bond_w.Uangle, host.E_bond_w.Utor, host.E_bond_w.Uimp);
    printf("qp\t%f\t%f\t%f\t%f\n", host.E_bond_q.Ubond, host.E_bond_q.Uangle, host.E_bond_q.Utor, host.E_bond_q.Uimp);
    printf("\n");

    printf("[nonbonded]\n");
    printf("type\tUcoul\tUvdw\n");
    printf("pp\t%f\t%f\n", host.E_nonbond_pp.Ucoul, host.E_nonbond_pp.Uvdw);
    printf("pw\t%f\t%f\n", host.E_nonbond_pw.Ucoul, host.E_nonbond_pw.Uvdw);
    printf("ww\t%f\t%f\n", host.E_nonbond_ww.Ucoul, host.E_nonbond_ww.Uvdw);
    printf("qx\t%f\t%f\n", host.E_nonbond_qx.Ucoul, host.E_nonbond_qx.Uvdw);
    printf("\n");

    printf("[restraint]\n");
    printf("Uradx\t%f\n", host.E_restraint.Uradx);
    printf("Upolx\t%f\n", host.E_restraint.Upolx);
    printf("Ushell\t%f\n", host.E_restraint.Ushell);
    printf("Ufix\t%f\n", host.E_restraint.Ufix);
    printf("Upres\t%f\n", host.E_restraint.Upres);
    printf("Total\t%f\n", host.E_restraint.Urestr);
    printf("\n");

    printf("[q-energies]\n");
    printf("lambda\tSUM\tUbond\tUangle\tUtor\tUimp\tUcoul\tUvdw\tUrestr\n");
    for (int state = 0; state < host.n_lambdas; state++) {
        printf("%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\n", lambdas[state], host.EQ_total[state].Utot, host.EQ_bond[state].Ubond,
               host.EQ_bond[state].Uangle, host.EQ_bond[state].Utor, host.EQ_bond[state].Uimp, host.EQ_nonbond_qx[state].Ucoul, host.EQ_nonbond_qx[state].Uvdw, host.EQ_restraint[state].Urestr);
    }
    printf("\n");

    printf("[total]\n");
    printf("Ukin\t%f\n", host.E_total.Ukin);
    printf("Upot\t%f\n", host.E_total.Upot);
    printf("Utot\t%f\n", host.E_total.Utot);
    printf("\n");
}

// Write header (number of atoms) to output file
void write_header(const char* filename) {
    auto& ctx = Context::instance();
    FILE* fp;

    char path[1024];
    sprintf(path, "%s/output/%s", ctx.base_folder.c_str(), filename);

    fp = fopen(path, "w");

    fprintf(fp, "%d\n", ctx.n_atoms);

    fclose(fp);
}

// Write header (number of atoms & lambdas) to output file
void write_energy_header() {
    auto& ctx = Context::instance();
    auto *lambdas = ctx.lambdas->cpu_data_p;
    FILE* fp;

    char path[1024];
    sprintf(path, "%s/output/%s", ctx.base_folder.c_str(), "energies.csv");

    fp = fopen(path, "w");

    fprintf(fp, "%d\n", ctx.n_atoms);

    fprintf(fp, "lambdas\n");
    fprintf(fp, "%d\n", ctx.n_lambdas);

    for (int state = 0; state < ctx.n_lambdas; state++) {
        fprintf(fp, "%f\n", lambdas[state]);
    }

    fclose(fp);
}

// Write step number, coordinates of atoms to coordinate output file
void write_coords(int iteration) {
    auto& ctx = Context::instance();
    auto &coords = ctx.coords->cpu_data_p;
    if (iteration % ctx.md.trajectory != 0) return;
    FILE* fp;
    int i;

    char path[1024];
    sprintf(path, "%s/output/%s", ctx.base_folder.c_str(), "coords.csv");

    fp = fopen(path, "a");

    fprintf(fp, "%d\n", iteration / ctx.md.trajectory);
    for (i = 0; i < ctx.n_atoms; i++) {
        fprintf(fp, "%f;%f;%f\n", coords[i].x, coords[i].y, coords[i].z);
    }

    fclose(fp);
}

// Write step number, velocities of atoms to coordinate output file
void write_velocities(int iteration) {
    auto& ctx = Context::instance();
    auto &velocities = ctx.velocities->cpu_data_p;
    if (iteration % ctx.md.trajectory != 0) return;
    FILE* fp;
    int i;

    char path[1024];
    sprintf(path, "%s/output/%s", ctx.base_folder.c_str(), "velocities.csv");

    fp = fopen(path, "a");

    fprintf(fp, "%d\n", iteration / ctx.md.trajectory);
    for (i = 0; i < ctx.n_atoms; i++) {
        fprintf(fp, "%f;%f;%f\n", velocities[i].x, velocities[i].y, velocities[i].z);
    }

    fclose(fp);
}

// Write step number, energies of atoms to coordinate output file
void write_energies(int iteration) {
    auto& ctx = Context::instance();
    auto *lambdas = ctx.lambdas->cpu_data_p;
    if (iteration % ctx.md.energy != 0) return;
    FILE* fp;

    char path[1024];
    sprintf(path, "%s/output/%s", ctx.base_folder.c_str(), "energies.csv");

    fp = fopen(path, "a");

    fprintf(fp, "interval %d\n", iteration / ctx.md.energy);

    fprintf(fp, "[temperature]\n");
    fprintf(fp, "Temp\t%f\n", ctx.Temp);
    fprintf(fp, "\n");

    fprintf(fp, "[bonded]\n");
    fprintf(fp, "type\tUbond\tUangle\tUtor\tUimp\n");
    fprintf(fp, "p\t%f\t%f\t%f\t%f\n", ctx.E_bond_p.Ubond, ctx.E_bond_p.Uangle, ctx.E_bond_p.Utor, ctx.E_bond_p.Uimp);
    fprintf(fp, "w\t%f\t%f\t%f\t%f\n", ctx.E_bond_w.Ubond, ctx.E_bond_w.Uangle, ctx.E_bond_w.Utor, ctx.E_bond_w.Uimp);
    fprintf(fp, "qp\t%f\t%f\t%f\t%f\n", ctx.E_bond_q.Ubond, ctx.E_bond_q.Uangle, ctx.E_bond_q.Utor, ctx.E_bond_q.Uimp);
    fprintf(fp, "\n");

    fprintf(fp, "[nonbonded]\n");
    fprintf(fp, "type\tUcoul\tUvdw\n");
    fprintf(fp, "pp\t%f\t%f\n", ctx.E_nonbond_pp.Ucoul, ctx.E_nonbond_pp.Uvdw);
    fprintf(fp, "pw\t%f\t%f\n", ctx.E_nonbond_pw.Ucoul, ctx.E_nonbond_pw.Uvdw);
    fprintf(fp, "ww\t%f\t%f\n", ctx.E_nonbond_ww.Ucoul, ctx.E_nonbond_ww.Uvdw);
    fprintf(fp, "qx\t%f\t%f\n", ctx.E_nonbond_qx.Ucoul, ctx.E_nonbond_qx.Uvdw);
    fprintf(fp, "\n");

    fprintf(fp, "[restraint]\n");
    fprintf(fp, "Uradx\t%f\n", ctx.E_restraint.Uradx);
    fprintf(fp, "Upolx\t%f\n", ctx.E_restraint.Upolx);
    fprintf(fp, "Ushell\t%f\n", ctx.E_restraint.Ushell);
    fprintf(fp, "Ufix\t%f\n", ctx.E_restraint.Ufix);
    fprintf(fp, "Upres\t%f\n", ctx.E_restraint.Upres);
    fprintf(fp, "Total\t%f\n", ctx.E_restraint.Urestr);
    fprintf(fp, "\n");

    fprintf(fp, "[q-energies]\n");
    fprintf(fp, "lambda\tSUM\tUbond\tUangle\tUtor\tUimp\tUcoul\tUvdw\tUrestr\n");
    for (int state = 0; state < ctx.n_lambdas; state++) {
        fprintf(fp, "%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\n", lambdas[state], ctx.EQ_total[state].Utot, ctx.EQ_bond[state].Ubond,
                ctx.EQ_bond[state].Uangle, ctx.EQ_bond[state].Utor, ctx.EQ_bond[state].Uimp, ctx.EQ_nonbond_qx[state].Ucoul, ctx.EQ_nonbond_qx[state].Uvdw, ctx.EQ_restraint[state].Urestr);
    }
    fprintf(fp, "\n");

    fprintf(fp, "[total]\n");
    fprintf(fp, "Ukin\t%f\n", ctx.E_total.Ukin);
    fprintf(fp, "Upot\t%f\n", ctx.E_total.Upot);
    fprintf(fp, "Utot\t%f\n", ctx.E_total.Utot);
    fprintf(fp, "\n");

    fclose(fp);
}
