#include "csv_out.h"

#include <cstdio>
#include <filesystem>
#include <stdexcept>
#include <string>

namespace {
FILE* open_file_or_throw(const std::string& path, const char* mode) {
    FILE* fp = std::fopen(path.c_str(), mode);
    if (!fp) {
        throw std::runtime_error("Could not open CSV output file " + path);
    }
    return fp;
}
}  // namespace

CsvOutput::CsvOutput(const std::string& csv_dir) : csv_dir_(csv_dir) {}

std::string CsvOutput::output_path(const char* filename) const {
    return (std::filesystem::path(csv_dir_) / "output" / filename).string();
}

void CsvOutput::write_atom_count_header(Context& ctx, const char* filename) const {
    FILE* fp = open_file_or_throw(output_path(filename), "w");
    std::fprintf(fp, "%d\n", ctx.n_atoms);
    std::fclose(fp);
}

void CsvOutput::write_energy_header(Context& ctx) const {
    FILE* fp = open_file_or_throw(output_path("energies.csv"), "w");
    auto* lambdas = ctx.lambdas ? ctx.lambdas->cpu_data_p : nullptr;

    std::fprintf(fp, "%d\n", ctx.n_atoms);
    std::fprintf(fp, "lambdas\n");
    std::fprintf(fp, "%d\n", ctx.n_lambdas);
    for (int state = 0; state < ctx.n_lambdas; state++) {
        std::fprintf(fp, "%f\n", lambdas[state]);
    }

    std::fclose(fp);
}

void CsvOutput::init(Context& ctx) {
    std::filesystem::create_directories(std::filesystem::path(csv_dir_) / "output");
    write_atom_count_header(ctx, "coords.csv");
    write_atom_count_header(ctx, "velocities.csv");
    write_energy_header(ctx);
}

void CsvOutput::output_trajectory(Context& ctx, int iteration) {
    if (ctx.md.trajectory <= 0) return;
    if (iteration % ctx.md.trajectory != 0) return;

    FILE* coords_fp = open_file_or_throw(output_path("coords.csv"), "a");
    std::fprintf(coords_fp, "%d\n", iteration / ctx.md.trajectory);
    auto& coords = ctx.coords->cpu_data_p;
    for (int i = 0; i < ctx.n_atoms; i++) {
        std::fprintf(coords_fp, "%f;%f;%f\n", coords[i].x, coords[i].y, coords[i].z);
    }
    std::fclose(coords_fp);

    FILE* velocities_fp = open_file_or_throw(output_path("velocities.csv"), "a");
    std::fprintf(velocities_fp, "%d\n", iteration / ctx.md.trajectory);
    auto& velocities = ctx.velocities->cpu_data_p;
    for (int i = 0; i < ctx.n_atoms; i++) {
        std::fprintf(velocities_fp, "%f;%f;%f\n", velocities[i].x, velocities[i].y, velocities[i].z);
    }
    std::fclose(velocities_fp);
}

void CsvOutput::output_energy(Context& ctx, int iteration) {
    if (ctx.md.energy <= 0) return;
    if (iteration % ctx.md.energy != 0) return;

    FILE* fp = open_file_or_throw(output_path("energies.csv"), "a");
    auto* lambdas = ctx.lambdas ? ctx.lambdas->cpu_data_p : nullptr;
    auto* EQ_restraint = ctx.EQ_restraint->cpu_data_p;

    std::fprintf(fp, "interval %d\n", iteration / ctx.md.energy);

    std::fprintf(fp, "[temperature]\n");
    std::fprintf(fp, "Temp\t%f\n", ctx.Temp);
    std::fprintf(fp, "\n");

    std::fprintf(fp, "[bonded]\n");
    std::fprintf(fp, "type\tUbond\tUangle\tUtor\tUimp\n");
    std::fprintf(fp, "p\t%f\t%f\t%f\t%f\n", ctx.E_bond_p.Ubond, ctx.E_bond_p.Uangle, ctx.E_bond_p.Utor, ctx.E_bond_p.Uimp);
    std::fprintf(fp, "w\t%f\t%f\t%f\t%f\n", ctx.E_bond_w.Ubond, ctx.E_bond_w.Uangle, ctx.E_bond_w.Utor, ctx.E_bond_w.Uimp);
    std::fprintf(fp, "qp\t%f\t%f\t%f\t%f\n", ctx.E_bond_q.Ubond, ctx.E_bond_q.Uangle, ctx.E_bond_q.Utor, ctx.E_bond_q.Uimp);
    std::fprintf(fp, "\n");

    std::fprintf(fp, "[nonbonded]\n");
    std::fprintf(fp, "type\tUcoul\tUvdw\n");
    std::fprintf(fp, "pp\t%f\t%f\n", ctx.E_nonbond_pp.Ucoul, ctx.E_nonbond_pp.Uvdw);
    std::fprintf(fp, "pw\t%f\t%f\n", ctx.E_nonbond_pw.Ucoul, ctx.E_nonbond_pw.Uvdw);
    std::fprintf(fp, "ww\t%f\t%f\n", ctx.E_nonbond_ww.Ucoul, ctx.E_nonbond_ww.Uvdw);
    std::fprintf(fp, "qx\t%f\t%f\n", ctx.E_nonbond_qx.Ucoul, ctx.E_nonbond_qx.Uvdw);
    std::fprintf(fp, "\n");

    std::fprintf(fp, "[restraint]\n");
    std::fprintf(fp, "Uradx\t%f\n", ctx.E_restraint.Uradx);
    std::fprintf(fp, "Upolx\t%f\n", ctx.E_restraint.Upolx);
    std::fprintf(fp, "Ushell\t%f\n", ctx.E_restraint.Ushell);
    std::fprintf(fp, "Ufix\t%f\n", ctx.E_restraint.Ufix);
    std::fprintf(fp, "Upres\t%f\n", ctx.E_restraint.Upres);
    std::fprintf(fp, "Total\t%f\n", ctx.E_restraint.Urestr);
    std::fprintf(fp, "\n");

    std::fprintf(fp, "[q-energies]\n");
    std::fprintf(fp, "lambda\tSUM\tUbond\tUangle\tUtor\tUimp\tUcoul\tUvdw\tUrestr\n");
    for (int state = 0; state < ctx.n_lambdas; state++) {
        std::fprintf(fp,
                     "%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\n",
                     lambdas[state],
                     ctx.EQ_total[state].Utot,
                     ctx.EQ_bond[state].Ubond,
                     ctx.EQ_bond[state].Uangle,
                     ctx.EQ_bond[state].Utor,
                     ctx.EQ_bond[state].Uimp,
                     ctx.EQ_nonbond_qx[state].Ucoul,
                     ctx.EQ_nonbond_qx[state].Uvdw,
                     EQ_restraint[state].Urestr);
    }
    std::fprintf(fp, "\n");

    std::fprintf(fp, "[total]\n");
    std::fprintf(fp, "Ukin\t%f\n", ctx.E_total.Ukin);
    std::fprintf(fp, "Upot\t%f\n", ctx.E_total.Upot);
    std::fprintf(fp, "Utot\t%f\n", ctx.E_total.Utot);
    std::fprintf(fp, "\n");

    std::fclose(fp);
}
