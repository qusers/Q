#include "csv_parser.h"

#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <filesystem>
#include <fstream>
#include <stdexcept>
#include <string>
#include <vector>

#include "constants.h"

namespace {
struct CsvRows {
    int count = 0;
    std::vector<std::vector<std::string>> rows;
};

const std::string& field_or_empty(const std::vector<std::string>& row, size_t index) {
    static const std::string empty;
    return index < row.size() ? row[index] : empty;
}

int row_int(const std::vector<std::string>& row, size_t index) {
    return to_int(field_or_empty(row, index));
}

double row_double(const std::vector<std::string>& row, size_t index) {
    return to_double(field_or_empty(row, index));
}

bool row_on(const std::vector<std::string>& row, size_t index) {
    return field_or_empty(row, index) == "on";
}

CsvRows read_csv_rows(const std::string& path, bool required = true) {
    std::ifstream input(path.c_str());
    if (!input) {
        if (!required) {
            return {};
        }
        throw std::runtime_error("Could not open CSV file " + path);
    }

    CsvRows csv;
    std::string line;
    if (!std::getline(input, line)) {
        return csv;
    }
    csv.count = to_int(strip_cr(line));

    while (std::getline(input, line)) {
        csv.rows.push_back(split_semicolon(line));
    }

    return csv;
}

coord_t coord_from_row(const std::vector<std::string>& row) {
    return {row_double(row, 0), row_double(row, 1), row_double(row, 2)};
}

vel_t vel_from_row(const std::vector<std::string>& row) {
    return {row_double(row, 0), row_double(row, 1), row_double(row, 2)};
}

int lambda_count(const ParseResult& result) {
    return static_cast<int>(result.lambdas.size());
}
}  // namespace

CsvParser::CsvParser(const std::string& csv_dir) : BaseParser(csv_dir) {}

std::string CsvParser::path_for(const std::string& filename) const {
    std::filesystem::path base(file_path);
    if (std::filesystem::is_regular_file(base)) {
        base = base.parent_path();
    }
    return (base / filename).string();
}

void CsvParser::parse_md() {
    CsvRows file = read_csv_rows(path_for("md.csv"));
    md_t& md = result.md;

    md.steps = row_int(file.rows[0], 1);
    md.stepsize = row_double(file.rows[1], 1);
    md.temperature = row_double(file.rows[2], 1);

    md.thermostat = field_or_empty(file.rows[3], 1);

    md.bath_coupling = row_double(file.rows[4], 1);
    md.random_seed = row_int(file.rows[5], 1);
    md.initial_temperature = row_double(file.rows[6], 1);
    md.shake_solvent = row_on(file.rows[7], 1);
    md.shake_solute = row_on(file.rows[8], 1);
    md.shake_hydrogens = row_on(file.rows[9], 1);
    md.lrf = row_on(file.rows[10], 1);
    md.charge_groups = row_on(file.rows[11], 1);

    md.solute_solute = row_double(file.rows[12], 1);
    md.solvent_solvent = row_double(file.rows[13], 1);
    md.solute_solvent = row_double(file.rows[14], 1);
    md.q_atom = row_double(file.rows[15], 1);

    md.shell_radius = row_double(file.rows[16], 1);
    md.shell_force = row_double(file.rows[17], 1);

    md.radial_force = row_double(file.rows[18], 1);
    md.polarisation = true;
    md.charge_correction = md.polarisation;
    md.polarisation_force = row_double(file.rows[20], 1);

    md.non_bond = row_int(file.rows[21], 1);
    md.output = row_int(file.rows[22], 1);
    md.energy = row_int(file.rows[23], 1);
    md.trajectory = row_int(file.rows[24], 1);
}

void CsvParser::parse_lambdas() {
    CsvRows file = read_csv_rows(path_for("md.csv"));
    if (file.rows.size() <= 25) {
        return;
    }

    int index = 25;
    int n_lambdas = row_int(file.rows[index], 0);
    index++;
    result.lambdas.resize(n_lambdas);
    for (int i = 0; i < n_lambdas; i++, index++) {
        result.lambdas[i] = row_double(file.rows[index], 0);
    }

    int n_restrseqs = row_int(file.rows[index], 0);
    index++;
    result.restrseqs.resize(n_restrseqs);
    for (int i = 0; i < n_restrseqs; i++, index++) {
        result.restrseqs[i] = {
            row_int(file.rows[index], 0),
            row_int(file.rows[index], 1),
            row_double(file.rows[index], 2),
            row_int(file.rows[index], 3) == 1,
            row_int(file.rows[index], 4),
        };
    }

    int n_restrspos = row_int(file.rows[index], 0);
    index++;
    result.restrspos.resize(n_restrspos);
    for (int i = 0; i < n_restrspos; i++, index++) {
        result.restrspos[i] = {
            row_int(file.rows[index], 0),
            row_int(file.rows[index], 1),
            {row_double(file.rows[index], 2), row_double(file.rows[index], 3), row_double(file.rows[index], 4)},
            {row_double(file.rows[index], 5), row_double(file.rows[index], 6), row_double(file.rows[index], 7)},
        };
    }

    int n_restrdists = row_int(file.rows[index], 0);
    index++;
    result.restrdists.resize(n_restrdists);
    for (int i = 0; i < n_restrdists; i++, index++) {
        restrdis_t restrdist = {};
        restrdist.ai = row_int(file.rows[index], 0);
        restrdist.aj = row_int(file.rows[index], 1);
        restrdist.d1 = row_double(file.rows[index], 2);
        restrdist.d2 = row_double(file.rows[index], 3);
        restrdist.k = row_double(file.rows[index], 4);
        restrdist.ipsi = row_int(file.rows[index], 5);
        result.restrdists[i] = restrdist;
    }

    int n_restrangs = row_int(file.rows[index], 0);
    index++;
    result.restrangs.resize(n_restrangs);
    for (int i = 0; i < n_restrangs; i++, index++) {
        result.restrangs[i] = {
            row_int(file.rows[index], 0),
            row_int(file.rows[index], 1),
            row_int(file.rows[index], 2),
            row_int(file.rows[index], 3),
            row_double(file.rows[index], 4),
            row_double(file.rows[index], 5),
        };
    }

    int n_restrwalls = row_int(file.rows[index], 0);
    index++;
    result.restrwalls.resize(n_restrwalls);
    for (int i = 0; i < n_restrwalls; i++, index++) {
        result.restrwalls[i] = {
            row_int(file.rows[index], 0),
            row_int(file.rows[index], 1),
            row_double(file.rows[index], 2),
            row_double(file.rows[index], 3),
            row_double(file.rows[index], 5),
            row_double(file.rows[index], 4),
            row_int(file.rows[index], 6) == 1,
        };
    }
}

void CsvParser::parse_topo() {
    CsvRows file = read_csv_rows(path_for("topo.csv"));
    if (file.rows.empty()) {
        return;
    }

    result.topo.solvent_type = row_int(file.rows[0], 0);
    result.topo.exclusion_radius = row_double(file.rows[1], 0);
    result.topo.solvent_radius = row_double(file.rows[2], 0);
    result.topo.solute_center = coord_from_row(file.rows[3]);
    result.topo.solvent_center = coord_from_row(file.rows[4]);
    result.topo.el14_scale = row_double(file.rows[5], 0);
    result.topo.coulomb_constant = row_double(file.rows[6], 0);
    result.topo.vdw_rule = row_int(file.rows[7], 0);
    if (result.topo.vdw_rule == 0) {
        result.topo.vdw_rule = VDW_GEOMETRIC;
    }
    if (result.topo.vdw_rule != VDW_GEOMETRIC && result.topo.vdw_rule != VDW_ARITHMETIC) {
        throw std::runtime_error("Invalid vdw_rule in topo.csv");
    }
}

void CsvParser::parse_coords_init() {
    CsvRows file = read_csv_rows(path_for("coords.csv"));
    result.n_atoms_solute = file.rows.empty() ? 0 : row_int(file.rows[0], 0);
    result.coords_init.resize(file.count);

    for (int i = 0; i < file.count; i++) {
        result.coords_init[i] = coord_from_row(file.rows[i + 1]);
    }
}

void CsvParser::parse_coords() {
    const std::string coords_path = path_for("i_coords.csv");
    const std::string velocities_path = path_for("i_velocities.csv");
    const bool has_restart_coords = std::filesystem::exists(coords_path);
    const bool has_restart_velocities = std::filesystem::exists(velocities_path);
    if (has_restart_coords != has_restart_velocities) {
        throw std::runtime_error("CSV restart requires both i_coords.csv and i_velocities.csv.");
    }
    result.fresh_start = !has_restart_coords;

    CsvRows file = read_csv_rows(coords_path, false);
    if (file.count == 0) {
        result.coords = result.coords_init;
        return;
    }

    result.coords.resize(file.count);
    for (int i = 0; i < file.count; i++) {
        result.coords[i] = coord_from_row(file.rows[i]);
    }
}

void CsvParser::parse_velocities() {
    CsvRows file = read_csv_rows(path_for("i_velocities.csv"), false);
    if (file.count == 0) {
        return;
    }

    result.velocities.resize(file.count);
    for (int i = 0; i < file.count; i++) {
        result.velocities[i] = vel_from_row(file.rows[i]);
    }
}

void CsvParser::parse_bonds() {
    CsvRows file = read_csv_rows(path_for("bonds.csv"));
    result.n_bonds_solute = file.rows.empty() ? 0 : row_int(file.rows[0], 0);
    result.bonds.resize(file.count);
    for (int i = 0; i < file.count; i++) {
        result.bonds[i] = {row_int(file.rows[i + 1], 0), row_int(file.rows[i + 1], 1), row_int(file.rows[i + 1], 2)};
    }
}

void CsvParser::parse_cbonds() {
    CsvRows file = read_csv_rows(path_for("cbonds.csv"));
    result.cbonds.resize(file.count);
    for (int i = 0; i < file.count; i++) {
        result.cbonds[i] = {row_int(file.rows[i], 0), row_double(file.rows[i], 1), row_double(file.rows[i], 2)};
    }
}

void CsvParser::parse_angles() {
    CsvRows file = read_csv_rows(path_for("angles.csv"));
    result.n_angles_solute = file.rows.empty() ? 0 : row_int(file.rows[0], 0);
    result.angles.resize(file.count);
    for (int i = 0; i < file.count; i++) {
        result.angles[i] = {row_int(file.rows[i + 1], 0), row_int(file.rows[i + 1], 1), row_int(file.rows[i + 1], 2), row_int(file.rows[i + 1], 3)};
    }
}

void CsvParser::parse_cangles() {
    CsvRows file = read_csv_rows(path_for("cangles.csv"));
    result.cangles.resize(file.count);
    for (int i = 0; i < file.count; i++) {
        result.cangles[i] = {row_int(file.rows[i], 0), row_double(file.rows[i], 1), row_double(file.rows[i], 2)};
    }
}

void CsvParser::parse_torsions() {
    CsvRows file = read_csv_rows(path_for("torsions.csv"));
    result.n_torsions_solute = file.rows.empty() ? 0 : row_int(file.rows[0], 0);
    result.torsions.resize(file.count);
    for (int i = 0; i < file.count; i++) {
        result.torsions[i] = {row_int(file.rows[i + 1], 0), row_int(file.rows[i + 1], 1), row_int(file.rows[i + 1], 2), row_int(file.rows[i + 1], 3), row_int(file.rows[i + 1], 4)};
    }
}

void CsvParser::parse_ctorsions() {
    CsvRows file = read_csv_rows(path_for("ctorsions.csv"));
    result.ctorsions.resize(file.count);
    for (int i = 0; i < file.count; i++) {
        double paths = row_double(file.rows[i], 4);
        result.ctorsions[i] = {row_int(file.rows[i], 0), row_double(file.rows[i], 1), row_double(file.rows[i], 2), row_double(file.rows[i], 3), paths == 0.0 ? 0.0 : 1.0 / paths};
    }
}

void CsvParser::parse_impropers() {
    CsvRows file = read_csv_rows(path_for("impropers.csv"));
    result.n_impropers_solute = file.rows.empty() ? 0 : row_int(file.rows[0], 0);
    result.impropers.resize(file.count);
    for (int i = 0; i < file.count; i++) {
        result.impropers[i] = {row_int(file.rows[i + 1], 0), row_int(file.rows[i + 1], 1), row_int(file.rows[i + 1], 2), row_int(file.rows[i + 1], 3), row_int(file.rows[i + 1], 4)};
    }
}

void CsvParser::parse_cimpropers() {
    CsvRows file = read_csv_rows(path_for("cimpropers.csv"));
    result.cimpropers.resize(file.count);
    for (int i = 0; i < file.count; i++) {
        result.cimpropers[i] = {row_int(file.rows[i], 0), row_double(file.rows[i], 1), row_double(file.rows[i], 2)};
    }
}

void CsvParser::parse_restrspos() {}
void CsvParser::parse_restrangs() {}
void CsvParser::parse_restrdists() {}
void CsvParser::parse_restrseqs() {}
void CsvParser::parse_restrwalls() {}

void CsvParser::parse_charges() {
    CsvRows file = read_csv_rows(path_for("charges.csv"));
    result.charges.resize(file.count);
    for (int i = 0; i < file.count; i++) {
        result.charges[i] = {row_int(file.rows[i], 0), row_int(file.rows[i], 1)};
    }
}

void CsvParser::parse_ccharges() {
    CsvRows file = read_csv_rows(path_for("ccharges.csv"));
    result.ccharges.resize(file.count);
    for (int i = 0; i < file.count; i++) {
        result.ccharges[i] = {row_int(file.rows[i], 0), row_double(file.rows[i], 1)};
    }
}

void CsvParser::parse_atypes() {
    CsvRows file = read_csv_rows(path_for("atypes.csv"));
    result.atypes.resize(file.count);
    for (int i = 0; i < file.count; i++) {
        result.atypes[i] = {row_int(file.rows[i], 0), row_int(file.rows[i], 1)};
    }
}

void CsvParser::parse_catypes() {
    CsvRows file = read_csv_rows(path_for("catypes.csv"));
    result.catypes.resize(file.count);
    for (int i = 0; i < file.count; i++) {
        result.catypes[i] = {
            row_int(file.rows[i], 0),
            row_double(file.rows[i], 1),
            row_double(file.rows[i], 2),
            row_double(file.rows[i], 3),
            row_double(file.rows[i], 6),
            row_double(file.rows[i], 7),
        };
    }
}

void CsvParser::parse_heavy() {}

void CsvParser::parse_excluded() {
    std::ifstream input(path_for("excluded.csv").c_str());
    if (!input) {
        throw std::runtime_error("Could not open CSV file " + path_for("excluded.csv"));
    }

    std::string line;
    std::getline(input, line);
    result.excluded.assign(result.coords_init.size(), false);
    for (size_t i = 0; i < result.excluded.size() && i < line.size(); i++) {
        result.excluded[i] = line[i] == '1';
    }
}

void CsvParser::parse_ngbrs14() {
    CsvRows file = read_csv_rows(path_for("ngbrs14.csv"));
    result.ngbrs14.clear();
    int n_lines = file.count;
    for (int line_i = 0; line_i < n_lines && line_i < static_cast<int>(file.rows.size()); line_i++) {
        const std::string& line = field_or_empty(file.rows[line_i], 0);
        for (int i = 0; i < line_width && i < static_cast<int>(line.size()); i++) {
            if (line[i] == '1') {
                result.ngbrs14.push_back({line_i, (line_i + i + 1) % n_lines});
            }
        }
    }
}

void CsvParser::parse_ngbrs14_long() {
    CsvRows file = read_csv_rows(path_for("ngbrs14long.csv"));
    result.ngbrs14_long.resize(file.count);
    for (int i = 0; i < file.count; i++) {
        result.ngbrs14_long[i] = {row_int(file.rows[i], 0) - 1, row_int(file.rows[i], 1) - 1};
    }
}

void CsvParser::parse_ngbrs23() {
    CsvRows file = read_csv_rows(path_for("ngbrs23.csv"));
    result.ngbrs23.clear();
    int n_lines = file.count;
    for (int line_i = 0; line_i < n_lines && line_i < static_cast<int>(file.rows.size()); line_i++) {
        const std::string& line = field_or_empty(file.rows[line_i], 0);
        for (int i = 0; i < line_width && i < static_cast<int>(line.size()); i++) {
            if (line[i] == '1') {
                result.ngbrs23.push_back({line_i, (line_i + i + 1) % n_lines});
            }
        }
    }
}

void CsvParser::parse_ngbrs23_long() {
    CsvRows file = read_csv_rows(path_for("ngbrs23long.csv"));
    result.ngbrs23_long.resize(file.count);
    for (int i = 0; i < file.count; i++) {
        result.ngbrs23_long[i] = {row_int(file.rows[i], 0) - 1, row_int(file.rows[i], 1) - 1};
    }
}

void CsvParser::parse_molecules() {
    CsvRows file = read_csv_rows(path_for("molecules.csv"));
    result.molecules.resize(file.count);
    for (int i = 0; i < file.count; i++) {
        result.molecules[i] = row_int(file.rows[i], 0);
    }
}

void CsvParser::parse_charge_groups() {
    if (!result.md.charge_groups) {
        return;
    }

    CsvRows file = read_csv_rows(path_for("charge_groups.csv"));
    if (file.rows.empty()) {
        return;
    }

    charge_group_config_t config;
    config.n_cgrps_solute = row_int(file.rows[0], 0);
    config.n_cgrps_solvent = row_int(file.rows[0], 1);
    config.iuse_switch_atom = row_int(file.rows[0], 2);

    int n_charge_groups = config.n_cgrps_solute + config.n_cgrps_solvent;
    config.charge_groups.resize(n_charge_groups);

    int row = 1;
    for (int i = 0; i < n_charge_groups; i++) {
        int n_atoms = row_int(file.rows[row], 0);
        config.charge_groups[i].iswitch = row_int(file.rows[row], 1);
        config.charge_groups[i].atoms.resize(n_atoms);
        row++;
        for (int j = 0; j < n_atoms; j++, row++) {
            config.charge_groups[i].atoms[j] = row_int(file.rows[row], 0);
        }
    }

    result.charge_groups.clear();
    result.charge_groups.push_back(config);
}

void CsvParser::parse_q_atoms() {
    CsvRows file = read_csv_rows(path_for("q_atoms.csv"));
    result.q_atoms.resize(file.count);
    for (int i = 0; i < file.count; i++) {
        result.q_atoms[i] = row_int(file.rows[i], 0) - 1;
    }
}

void CsvParser::parse_q_angcouples() {
    CsvRows file = read_csv_rows(path_for("q_angcouples.csv"));
    result.q_angcouples.resize(file.count);
    for (int i = 0; i < file.count; i++) {
        result.q_angcouples[i] = {row_int(file.rows[i], 0), row_int(file.rows[i], 1)};
    }
}

void CsvParser::parse_q_atypes() {
    CsvRows file = read_csv_rows(path_for("q_atypes.csv"));
    int n_lambdas = lambda_count(result);
    result.q_atypes.resize(result.q_atoms.size() * n_lambdas);
    for (size_t i = 0; i < result.q_atoms.size(); i++) {
        for (int state = 0; state < n_lambdas; state++) {
            result.q_atypes[i + state * result.q_atoms.size()].code = row_int(file.rows[i + state * result.q_atoms.size()], 0);
        }
    }
}

void CsvParser::parse_q_catypes() {
    CsvRows file = read_csv_rows(path_for("q_catypes.csv"));
    result.q_catypes.resize(file.count);
    for (int i = 0; i < file.count; i++) {
        result.q_catypes[i] = {
            i + 1,
            row_double(file.rows[i], 7),
            row_double(file.rows[i], 1),
            row_double(file.rows[i], 2),
            row_double(file.rows[i], 5),
            row_double(file.rows[i], 6),
        };
    }
}

void CsvParser::parse_q_charges() {
    CsvRows file = read_csv_rows(path_for("q_charges.csv"));
    int n_lambdas = lambda_count(result);
    result.q_charges.resize(result.q_atoms.size() * n_lambdas);
    for (size_t i = 0; i < result.q_atoms.size(); i++) {
        for (int state = 0; state < n_lambdas; state++) {
            result.q_charges[i + state * result.q_atoms.size()].charge = row_double(file.rows[i + state * result.q_atoms.size()], 0);
        }
    }
}

void CsvParser::parse_q_angles() {
    CsvRows file = read_csv_rows(path_for("q_angles.csv"));
    int n_lambdas = lambda_count(result);
    int n_qangles = n_lambdas > 0 ? file.count / n_lambdas : 0;
    result.q_angles.resize(n_qangles * n_lambdas);
    for (int i = 0; i < n_qangles; i++) {
        for (int state = 0; state < n_lambdas; state++) {
            const auto& row = file.rows[i + state * n_qangles];
            result.q_angles[i + state * n_qangles] = {row_int(row, 0), row_int(row, 1), row_int(row, 2), row_int(row, 3)};
        }
    }
}

void CsvParser::parse_q_cangles() {
    CsvRows file = read_csv_rows(path_for("q_cangles.csv"));
    result.q_cangles.resize(file.count);
    for (int i = 0; i < file.count; i++) {
        result.q_cangles[i] = {0, row_double(file.rows[i], 0), row_double(file.rows[i], 1)};
    }
}

void CsvParser::parse_q_bonds() {
    CsvRows file = read_csv_rows(path_for("q_bonds.csv"));
    int n_lambdas = lambda_count(result);
    int n_qbonds = n_lambdas > 0 ? file.count / n_lambdas : 0;
    result.q_bonds.resize(n_qbonds * n_lambdas);
    for (int i = 0; i < n_qbonds; i++) {
        for (int state = 0; state < n_lambdas; state++) {
            const auto& row = file.rows[i + state * n_qbonds];
            result.q_bonds[i + state * n_qbonds] = {row_int(row, 0), row_int(row, 1), row_int(row, 2)};
        }
    }
}

void CsvParser::parse_q_cbonds() {
    CsvRows file = read_csv_rows(path_for("q_cbonds.csv"));
    result.q_cbonds.resize(file.count);
    for (int i = 0; i < file.count; i++) {
        result.q_cbonds[i] = {0, row_double(file.rows[i], 0), row_double(file.rows[i], 1)};
    }
}

void CsvParser::parse_q_impropers() {
    CsvRows file = read_csv_rows(path_for("q_impropers.csv"));
    int n_lambdas = lambda_count(result);
    int n_qimpropers = n_lambdas > 0 ? file.count / n_lambdas : 0;
    result.q_impropers.resize(n_qimpropers * n_lambdas);
    for (int i = 0; i < n_qimpropers; i++) {
        for (int state = 0; state < n_lambdas; state++) {
            const auto& row = file.rows[i + state * n_qimpropers];
            result.q_impropers[i + state * n_qimpropers] = {row_int(row, 0), row_int(row, 1), row_int(row, 2), row_int(row, 3), row_int(row, 4)};
        }
    }
}

void CsvParser::parse_q_cimpropers() {
    CsvRows file = read_csv_rows(path_for("q_cimpropers.csv"));
    result.q_cimpropers.resize(file.count);
    for (int i = 0; i < file.count; i++) {
        result.q_cimpropers[i] = {row_double(file.rows[i], 0), row_double(file.rows[i], 1)};
    }
}

void CsvParser::parse_q_torsions() {
    CsvRows file = read_csv_rows(path_for("q_torsions.csv"));
    int n_lambdas = lambda_count(result);
    int n_qtorsions = n_lambdas > 0 ? file.count / n_lambdas : 0;
    result.q_torsions.resize(n_qtorsions * n_lambdas);
    for (int i = 0; i < n_qtorsions; i++) {
        for (int state = 0; state < n_lambdas; state++) {
            const auto& row = file.rows[i + state * n_qtorsions];
            result.q_torsions[i + state * n_qtorsions] = {row_int(row, 0), row_int(row, 1), row_int(row, 2), row_int(row, 3), row_int(row, 4)};
        }
    }
}

void CsvParser::parse_q_ctorsions() {
    CsvRows file = read_csv_rows(path_for("q_ctorsions.csv"));
    result.q_ctorsions.resize(file.count);
    for (int i = 0; i < file.count; i++) {
        result.q_ctorsions[i] = {0, row_double(file.rows[i], 0), row_double(file.rows[i], 1), row_double(file.rows[i], 2), 0.0};
    }
}

void CsvParser::parse_q_offdiags() {
    CsvRows file = read_csv_rows(path_for("q_offdiags.csv"));
    result.q_offdiags.resize(file.count);
    for (int i = 0; i < file.count; i++) {
        result.q_offdiags[i] = {row_int(file.rows[i], 0), row_int(file.rows[i], 1), row_int(file.rows[i], 2), row_int(file.rows[i], 3), row_double(file.rows[i], 4), row_double(file.rows[i], 5)};
    }
}

void CsvParser::parse_q_imprcouples() {
    CsvRows file = read_csv_rows(path_for("q_imprcouples.csv"));
    result.q_imprcouples.resize(file.count);
    for (int i = 0; i < file.count; i++) {
        result.q_imprcouples[i] = {row_int(file.rows[i], 0), row_int(file.rows[i], 1)};
    }
}

void CsvParser::parse_q_softpairs() {
    CsvRows file = read_csv_rows(path_for("q_softpairs.csv"));
    result.q_softpairs.resize(file.count);
    for (int i = 0; i < file.count; i++) {
        result.q_softpairs[i] = {row_int(file.rows[i], 0), row_int(file.rows[i], 1)};
    }
}

void CsvParser::parse_q_torcouples() {
    CsvRows file = read_csv_rows(path_for("q_torcouples.csv"));
    result.q_torcouples.resize(file.count);
    for (int i = 0; i < file.count; i++) {
        result.q_torcouples[i] = {row_int(file.rows[i], 0), row_int(file.rows[i], 1)};
    }
}

void CsvParser::parse_q_elscales() {
    CsvRows file = read_csv_rows(path_for("q_elscales.csv"));
    int n_lambdas = lambda_count(result);
    int n_qelscales = n_lambdas > 0 ? file.count / n_lambdas : 0;
    result.q_elscales.resize(n_qelscales * n_lambdas);
    for (int i = 0; i < n_qelscales; i++) {
        for (int state = 0; state < n_lambdas; state++) {
            const auto& row = file.rows[i + state * n_qelscales];
            result.q_elscales[i + state * n_qelscales] = {row_int(row, 0), row_int(row, 1), row_double(row, 2)};
        }
    }
}

void CsvParser::parse_q_exclpairs() {
    CsvRows file = read_csv_rows(path_for("q_exclpairs.csv"));
    int n_lambdas = lambda_count(result);
    int n_qexclpairs = n_lambdas > 0 ? file.count / n_lambdas : 0;
    result.q_exclpairs.resize(n_qexclpairs * n_lambdas);
    for (int i = 0; i < n_qexclpairs; i++) {
        for (int state = 0; state < n_lambdas; state++) {
            const auto& row = file.rows[i + state * n_qexclpairs];
            result.q_exclpairs[i + state * n_qexclpairs] = {row_int(row, 0), row_int(row, 1), row_int(row, 2)};
        }
    }
}

void CsvParser::parse_q_shakes() {
    CsvRows file = read_csv_rows(path_for("q_shakes.csv"));
    int n_lambdas = lambda_count(result);
    int n_qshakes = n_lambdas > 0 ? file.count / n_lambdas : 0;
    result.q_shakes.resize(n_qshakes * n_lambdas);
    for (int i = 0; i < n_qshakes; i++) {
        for (int state = 0; state < n_lambdas; state++) {
            const auto& row = file.rows[i + state * n_qshakes];
            result.q_shakes[i + state * n_qshakes] = {row_int(row, 0), row_int(row, 1), row_double(row, 2)};
        }
    }
}

void CsvParser::parse_q_softcores() {
    CsvRows file = read_csv_rows(path_for("q_softcores.csv"));
    int n_lambdas = lambda_count(result);
    int n_qsoftcores = n_lambdas > 0 ? file.count / n_lambdas : 0;
    result.q_softcores.resize(n_qsoftcores * n_lambdas);
    for (int i = 0; i < n_qsoftcores; i++) {
        for (int state = 0; state < n_lambdas; state++) {
            result.q_softcores[i + state * n_qsoftcores].s = row_double(file.rows[i + state * n_qsoftcores], 0);
        }
    }
}
