#include "inp_parser.h"

#include <algorithm>
#include <cctype>
#include <cmath>
#include <cstdint>
#include <cstdlib>
#include <cstring>
#include <filesystem>
#include <fstream>
#include <map>
#include <stdexcept>
#include <string>
#include <vector>

#include "constants.h"
#include "str_helpers.h"

namespace {
std::runtime_error parse_error(const std::string& message) {
    return std::runtime_error(message);
}

std::string value_or(const std::map<std::string, std::string>& values, const std::string& key, const std::string& fallback) {
    auto it = values.find(key);
    if (it == values.end()) return fallback;
    return it->second;
}

std::string bool_value(const std::map<std::string, std::string>& values, const std::string& key, const std::string& fallback) {
    auto it = values.find(key);
    if (it == values.end()) return fallback;
    std::string value = lower_normalized(it->second);
    if (value == "true" || value == "yes" || value == "1") return "on";
    if (value == "false" || value == "no" || value == "0") return "off";
    return value;
}

int parse_int(const std::string& value) {
    return std::atoi(value.c_str());
}

double parse_double(const std::string& value) {
    return std::atof(value.c_str());
}

int row_int(const std::vector<std::string>& row, size_t index, int fallback = 0) {
    return index < row.size() ? parse_int(row[index]) : fallback;
}

double row_double(const std::vector<std::string>& row, size_t index, double fallback = 0.0) {
    return index < row.size() ? parse_double(row[index]) : fallback;
}

bool is_on_value(const std::map<std::string, std::string>& values, const std::string& key, const std::string& fallback) {
    return bool_value(values, key, fallback) == "on";
}

std::vector<std::vector<std::string>> checked_split_groups(const std::vector<std::string>& flat, size_t width) {
    if (width == 0) throw parse_error("Invalid grouped input width.");
    if (flat.size() % width != 0) {
        throw parse_error("Invalid grouped input length: " + std::to_string(flat.size()) + " values cannot be split into groups of " + std::to_string(width) + ".");
    }

    std::vector<std::vector<std::string>> groups;
    groups.reserve(flat.size() / width);
    for (size_t i = 0; i + width <= flat.size(); i += width) {
        groups.emplace_back(flat.begin() + i, flat.begin() + i + width);
    }
    return groups;
}

coord_t row_coord(const std::vector<std::string>& row) {
    return {row_double(row, 0), row_double(row, 1), row_double(row, 2)};
}

vel_t row_vel(const std::vector<std::string>& row) {
    return {row_double(row, 0), row_double(row, 1), row_double(row, 2)};
}

std::string dirname_of_local(const std::string& path) {
    std::filesystem::path p(path);
    std::filesystem::path parent = p.parent_path();
    if (parent.empty()) return ".";
    return parent.string();
}

std::string resolve_existing_or_cwd_local(const std::string& value, const std::string& input_dir) {
    if (value.empty()) return "";
    std::filesystem::path path(value);
    if (path.is_absolute()) return path.string();

    std::filesystem::path input_relative = std::filesystem::path(input_dir) / path;
    if (std::filesystem::exists(input_relative)) return input_relative.string();
    return (std::filesystem::current_path() / path).string();
}

std::string resolve_output_path_local(const std::string& value, const std::string& input_dir) {
    if (value.empty()) return "";
    std::filesystem::path path(value);
    if (path.is_absolute()) return path.string();
    return (std::filesystem::path(input_dir) / path).string();
}

std::string join_space(const std::vector<std::string>& fields) {
    std::string out;
    for (size_t i = 0; i < fields.size(); i++) {
        if (i > 0) out += " ";
        out += fields[i];
    }
    return out;
}

std::string trajectory_atoms_setting(const std::map<std::string, std::vector<std::vector<std::string>>>& lists) {
    auto it = lists.find("trajectory-atoms");
    if (it == lists.end() || it->second.empty()) return "";
    if (it->second.size() > 1) {
        throw parse_error("Native QGPU output supports only one trajectory_atoms mask line.");
    }
    return join_space(it->second[0]);
}

std::vector<std::string> split_state_row(const std::string& value) {
    return split_ws(value);
}

void add_short_neighbors(const std::vector<std::string>& lines, std::vector<std::pair<int, int>>& out) {
    out.clear();
    int n_lines = static_cast<int>(lines.size());
    for (int line_i = 0; line_i < n_lines; line_i++) {
        for (int i = 0; i < line_width && static_cast<size_t>(i) < lines[line_i].size(); i++) {
            if (lines[line_i][i] == '1') {
                out.push_back({line_i, (line_i + i + 1) % n_lines});
            }
        }
    }
}

void add_long_neighbors(const std::vector<std::vector<std::string>>& rows, std::vector<std::pair<int, int>>& out) {
    out.resize(rows.size());
    for (size_t i = 0; i < rows.size(); i++) {
        out[i] = {row_int(rows[i], 0) - 1, row_int(rows[i], 1) - 1};
    }
}

bool read_fortran_record(std::ifstream& in, std::vector<char>& payload) {
    int32_t nbytes = 0;
    in.read(reinterpret_cast<char*>(&nbytes), sizeof(nbytes));
    if (!in) return false;
    if (nbytes <= 0) throw parse_error("Invalid Fortran restart record marker.");
    payload.assign(nbytes, 0);
    in.read(payload.data(), nbytes);
    int32_t end_nbytes = 0;
    in.read(reinterpret_cast<char*>(&end_nbytes), sizeof(end_nbytes));
    if (!in || end_nbytes != nbytes) throw parse_error("Invalid or truncated Fortran restart record.");
    return true;
}

std::vector<std::vector<std::string>> unpack_restart_vector(const std::vector<char>& payload) {
    if (payload.size() < sizeof(int32_t)) throw parse_error("Invalid restart vector record.");
    int32_t nat3 = 0;
    std::memcpy(&nat3, payload.data(), sizeof(int32_t));
    if (nat3 <= 0 || nat3 % 3 != 0) throw parse_error("Invalid restart atom count.");
    size_t expected = sizeof(int32_t) + static_cast<size_t>(nat3) * sizeof(double);
    if (payload.size() != expected) throw parse_error("Invalid restart vector record length.");

    const double* values = reinterpret_cast<const double*>(payload.data() + sizeof(int32_t));
    std::vector<std::vector<std::string>> rows;
    rows.reserve(nat3 / 3);
    for (int i = 0; i < nat3; i += 3) {
        rows.push_back({std::to_string(values[i]), std::to_string(values[i + 1]), std::to_string(values[i + 2])});
    }
    return rows;
}

std::vector<real_t> unpack_restart_theta_corr(const std::vector<char>& payload) {
    if (payload.size() < sizeof(int32_t)) return {};

    int32_t n_shells = 0;
    std::memcpy(&n_shells, payload.data(), sizeof(int32_t));
    if (n_shells <= 0) return {};

    size_t expected = sizeof(int32_t) + static_cast<size_t>(n_shells) * sizeof(float);
    if (payload.size() != expected) return {};

    std::vector<real_t> theta_corr(n_shells);
    const float* values = reinterpret_cast<const float*>(payload.data() + sizeof(int32_t));
    for (int i = 0; i < n_shells; i++) {
        theta_corr[i] = static_cast<real_t>(values[i]);
    }
    return theta_corr;
}

void read_restart_vectors(const std::string& path,
                          std::vector<std::vector<std::string>>& coords,
                          std::vector<std::vector<std::string>>& velocities,
                          std::vector<real_t>& theta_corr) {
    std::ifstream in(path.c_str(), std::ios::binary);
    if (!in) throw parse_error("Could not open restart file " + path);

    std::vector<char> payload;
    if (!read_fortran_record(in, payload)) throw parse_error("Could not read restart coordinates.");
    coords = unpack_restart_vector(payload);
    if (!read_fortran_record(in, payload)) throw parse_error("Could not read restart velocities.");
    velocities = unpack_restart_vector(payload);
    if (coords.size() != velocities.size()) throw parse_error("Restart coordinate and velocity atom counts differ.");

    theta_corr.clear();
    if (read_fortran_record(in, payload)) {
        theta_corr = unpack_restart_theta_corr(payload);
    }
}

double q_randm(int& seed) {
    const int m = 100000000;
    const int m1 = 10000;
    const int mult = 31415821;
    int irand = std::abs(seed) % m;
    int irandh = irand / m1;
    int irandl = irand % m1;
    int multh = mult / m1;
    int multl = mult % m1;
    irand = ((irandh * multl + irandl * multh) % m1) * m1 + irandl * multl;
    irand = (irand + 1) % m;
    double r = ((irand / 10) * 10) / static_cast<double>(m);
    if (r <= 0.0 || r > 1.0) r = 0.0;
    seed = irand;
    return r;
}

double q_gauss(int& seed, double sd) {
    double total = 0.0;
    for (int i = 0; i < 12; i++) total += q_randm(seed);
    return (total - 6.0) * sd;
}
}  // namespace

struct InpParser::MDInput {
    std::map<std::string, std::map<std::string, std::string>> keyed;
    std::map<std::string, std::vector<std::vector<std::string>>> lists;
};

struct InpParser::TopData {
    std::string natoms_solute = "0";
    std::string nbonds_solute = "0";
    std::string nangles_solute = "0";
    std::string ntorsions_solute = "0";
    std::string nimpropers_solute = "0";
    std::vector<std::vector<std::string>> coords;
    std::vector<std::vector<std::string>> bonds;
    std::vector<std::vector<std::string>> angles;
    std::vector<std::vector<std::string>> torsions;
    std::vector<std::vector<std::string>> impropers;
    std::vector<std::vector<std::string>> ngbr14long;
    std::vector<std::vector<std::string>> ngbr23long;
    std::vector<std::pair<int, int>> atypes;
    std::map<int, std::vector<std::string>> cbonds;
    std::map<int, std::vector<std::string>> cangles;
    std::map<int, std::vector<std::string>> ctorsions;
    std::map<int, std::vector<std::string>> cimpropers;
    std::map<int, std::vector<std::string>> catypes;
    std::vector<std::pair<int, int>> charges;
    std::map<int, std::string> ccharges;
    std::vector<std::string> ngbr14;
    std::vector<std::string> ngbr23;
    std::vector<std::string> molecules;
    std::vector<std::string> excluded;
    std::string solute_cgps = "0";
    std::string solvent_cgps = "0";
    std::vector<std::vector<std::string>> charge_group_headers;
    std::vector<std::vector<std::string>> charge_group_atoms;
    std::string solvtype = "0";
    std::string exclusion = "0";
    std::string radii = "0";
    std::string scaling14 = "1";
    std::string coulomb = "332.0716";
    std::string vdw_rule = "1";
    std::vector<std::string> solucenter;
    std::vector<std::string> solvcenter;
};

struct InpParser::FepData {
    int states = 0;
    std::vector<std::string> q_atoms;
    std::vector<std::vector<std::string>> q_atypes;
    std::vector<std::vector<std::string>> q_charges;
    std::vector<std::vector<std::string>> q_softpairs;
    std::vector<std::vector<std::string>> q_exclpairs;
    std::vector<std::vector<std::string>> q_softcores;
    std::vector<std::vector<std::string>> q_bonds;
    std::vector<std::vector<std::string>> q_angles;
    std::vector<std::vector<std::string>> q_torsions;
    std::vector<std::vector<std::string>> q_impropers;
    std::vector<std::vector<std::string>> q_shakes;
    std::vector<std::vector<std::string>> q_catypes;
    std::vector<std::vector<std::string>> q_elscales;
    std::vector<std::vector<std::string>> q_cbonds;
    std::vector<std::vector<std::string>> q_cangles;
    std::vector<std::vector<std::string>> q_ctorsions;
    std::vector<std::vector<std::string>> q_cimpropers;
    std::vector<std::vector<std::string>> q_angcouples;
    std::vector<std::vector<std::string>> q_torcouples;
    std::vector<std::vector<std::string>> q_imprcouples;
    std::vector<std::vector<std::string>> q_offdiags;
    std::map<std::string, int> atypemap;
};

InpParser::InpParser(const std::string& inp_file) : BaseParser(inp_file), input_dir_(dirname_of_local(inp_file)) {}

InpParser::~InpParser() = default;

void InpParser::ensure_md_input() {
    if (input_) return;

    std::ifstream in(file_path.c_str());
    if (!in) throw parse_error("Could not open Q input file " + file_path);

    input_ = std::make_unique<MDInput>();
    std::string section;
    std::string raw;
    while (std::getline(in, raw)) {
        std::string line = strip_comment(raw);
        if (line.empty()) continue;
        if (line[0] == '[' && line[line.size() - 1] == ']') {
            section = lower_normalized(trim(line.substr(1, line.size() - 2)));
            continue;
        }
        if (section.empty()) continue;

        std::vector<std::string> fields = split_ws(line);
        if (fields.empty()) continue;

        if (section == "lambdas" || section == "trajectory-atoms" || section == "sequence-restraints" || section == "atom-restraints" ||
            section == "distance-restraints" || section == "angle-restraints" || section == "wall-restraints") {
            input_->lists[section].push_back(fields);
        } else if (fields.size() >= 2) {
            std::string key = lower_normalized(fields[0]);
            std::string value = fields[1];
            for (size_t i = 2; i < fields.size(); i++) value += " " + fields[i];
            input_->keyed[section][key] = value;
        }
    }

    if (!input_->keyed.count("files")) throw parse_error("Q input requires a [files] section.");
    const auto& files = input_->keyed["files"];
    topology_file_ = resolve_existing_or_cwd_local(value_or(files, "topology", ""), input_dir_);
    fep_file_ = resolve_existing_or_cwd_local(value_or(files, "fep", ""), input_dir_);
    restart_file_ = resolve_existing_or_cwd_local(value_or(files, "restart", ""), input_dir_);
    if (topology_file_.empty()) throw parse_error("Q input [files] section must specify topology.");
}

void InpParser::ensure_topology() {
    if (top_) return;
    ensure_md_input();

    std::ifstream in(topology_file_.c_str());
    if (!in) throw parse_error("Could not open topology file " + topology_file_);

    top_ = std::make_unique<TopData>();
    int block = 0;
    int atype_count = 0;
    int charge_group_switch = 0;
    int charge_group_atoms_expected = 0;
    std::vector<std::string> current_group_header;
    std::vector<std::string> current_group_atoms;
    std::vector<std::string> coord_flat, bond_flat, angle_flat, torsion_flat, improper_flat;
    std::vector<std::string> charges_tmp, masses, aii_normal, bii_normal, aii_polar, bii_polar, aii14, bii14;
    std::string ngbr14_flat, ngbr23_flat;

    std::string raw;
    while (std::getline(in, raw)) {
        std::string line = trim(raw);
        if (line.empty()) continue;

        if (line.find("Total no. of atoms") != std::string::npos) {
            std::vector<std::string> f = split_ws(line);
            if (f.size() > 1) top_->natoms_solute = f[1];
            block = 1;
            continue;
        }
        if (line.find("No. of integer atom codes") != std::string::npos) { block = 2; continue; }
        if (line.find("No. of bonds") != std::string::npos) {
            std::vector<std::string> f = split_ws(line);
            if (f.size() > 1) top_->nbonds_solute = f[1];
            block = 3;
            continue;
        }
        if (line.find("No. of bond codes") != std::string::npos) { block = 4; continue; }
        if (line.find("No. of angles") != std::string::npos) {
            std::vector<std::string> f = split_ws(line);
            if (f.size() > 1) top_->nangles_solute = f[1];
            block = 5;
            continue;
        }
        if (line.find("No. of angle codes") != std::string::npos) { block = 6; continue; }
        if (line.find("No. of torsions") != std::string::npos) {
            std::vector<std::string> f = split_ws(line);
            if (f.size() > 1) top_->ntorsions_solute = f[1];
            block = 7;
            continue;
        }
        if (line.find("No. of torsion codes") != std::string::npos) { block = 8; continue; }
        if (line.find("No. of impropers") != std::string::npos) {
            std::vector<std::string> f = split_ws(line);
            if (f.size() > 1) top_->nimpropers_solute = f[1];
            block = 9;
            continue;
        }
        if (line.find("No. of improper codes") != std::string::npos) { block = 10; continue; }
        if (line.find("No. of atomic charges") != std::string::npos) { block = 11; continue; }
        if (line.find("No. of charge groups") != std::string::npos) {
            std::vector<std::string> f = split_ws(line);
            int total = f.empty() ? 0 : parse_int(f[0]);
            int solute = f.size() > 1 ? parse_int(f[1]) : 0;
            top_->solute_cgps = std::to_string(solute);
            top_->solvent_cgps = std::to_string(total - solute);
            block = 12;
            charge_group_switch = 1;
            continue;
        }
        if (line.find("vdW combination rule") != std::string::npos) {
            std::vector<std::string> f = split_ws(line);
            if (!f.empty()) top_->vdw_rule = f[0];
            block = 13;
            continue;
        }
        if (line.find("Electrostatic 1-4 scaling factor") != std::string::npos) { block = 14; }
        if (line.find("Masses") != std::string::npos) { block = 15; continue; }
        if (line.find("sqrt (Aii) normal") != std::string::npos || line.find("R* normal:") != std::string::npos) { block = 16; continue; }
        if (line.find("sqrt (Bii) normal") != std::string::npos || line.find("epsilon normal:") != std::string::npos) { block = 17; continue; }
        if (line.find("sqrt (Aii) polar") != std::string::npos || line.find("R* polar:") != std::string::npos) { block = 18; continue; }
        if (line.find("sqrt (Bii) polar") != std::string::npos || line.find("epsilon polar:") != std::string::npos) { block = 19; continue; }
        if (line.find("sqrt (Aii) 1-4") != std::string::npos || line.find("R* 1-4:") != std::string::npos) { block = 20; continue; }
        if (line.find("sqrt (Bii) 1-4") != std::string::npos || line.find("epsilon 1-4:") != std::string::npos) { block = 21; continue; }
        if (line.find("No. of type-2 vdW interactions") != std::string::npos) { block = 22; continue; }
        if (line.find("No. of 1-4 neighbours") != std::string::npos) { block = 23; continue; }
        if (line.find("No. of long 1-4 nbrs") != std::string::npos) { block = 24; continue; }
        if (line.find("No. of exclusions") != std::string::npos) { block = 25; continue; }
        if (line.find("No. of long exclusions") != std::string::npos) { block = 26; continue; }
        if (line.find("No. of residues") != std::string::npos) { block = 27; continue; }
        if (line.find("Sequence") != std::string::npos) { block = 28; continue; }
        if (line.find("No. of separate molecules") != std::string::npos) { block = 29; continue; }
        if (line.find("No. of atom types") != std::string::npos) { block = 30; continue; }
        if (line.find("No. of SYBYL atom types") != std::string::npos) { block = 31; continue; }
        if (line.find("solvent type (0=SPC,1=3-atom,2=general)") != std::string::npos) {
            std::vector<std::string> f = split_ws(line);
            if (!f.empty()) top_->solvtype = f[0];
            if (top_->solvtype != "0") throw parse_error("solvent type other than SPC (0) is not supported");
            block = 32;
            continue;
        }
        if (line.find("No. of excluded atoms") != std::string::npos) { block = 33; continue; }

        std::vector<std::string> f = split_ws(line);
        switch (block) {
            case 1: coord_flat.insert(coord_flat.end(), f.begin(), f.end()); break;
            case 2:
                for (const std::string& value : f) top_->atypes.push_back({++atype_count, parse_int(value)});
                break;
            case 3: bond_flat.insert(bond_flat.end(), f.begin(), f.end()); break;
            case 4:
                if (f.size() >= 3) top_->cbonds[parse_int(f[0])] = std::vector<std::string>(f.begin() + 1, f.begin() + 3);
                break;
            case 5: angle_flat.insert(angle_flat.end(), f.begin(), f.end()); break;
            case 6:
                if (f.size() >= 3) top_->cangles[parse_int(f[0])] = std::vector<std::string>(f.begin() + 1, f.begin() + 3);
                break;
            case 7: torsion_flat.insert(torsion_flat.end(), f.begin(), f.end()); break;
            case 8:
                if (f.size() >= 5) top_->ctorsions[parse_int(f[0])] = std::vector<std::string>(f.begin() + 1, f.begin() + 5);
                break;
            case 9: improper_flat.insert(improper_flat.end(), f.begin(), f.end()); break;
            case 10:
                if (f.size() >= 3) top_->cimpropers[parse_int(f[0])] = std::vector<std::string>(f.begin() + 1, f.begin() + 3);
                break;
            case 11: charges_tmp.insert(charges_tmp.end(), f.begin(), f.end()); break;
            case 12:
                if (charge_group_switch == 1 && f.size() >= 2) {
                    current_group_header = f;
                    current_group_atoms.clear();
                    charge_group_atoms_expected = parse_int(f[0]);
                    charge_group_switch = 2;
                } else if (charge_group_switch == 2) {
                    current_group_atoms.insert(current_group_atoms.end(), f.begin(), f.end());
                    if (static_cast<int>(current_group_atoms.size()) >= charge_group_atoms_expected) {
                        top_->charge_group_headers.push_back(current_group_header);
                        top_->charge_group_atoms.push_back(current_group_atoms);
                        charge_group_switch = 1;
                    }
                }
                break;
            case 14:
                if (f.size() >= 2) {
                    top_->scaling14 = f[0];
                    top_->coulomb = f[1];
                }
                break;
            case 15: masses.insert(masses.end(), f.begin(), f.end()); break;
            case 16: aii_normal.insert(aii_normal.end(), f.begin(), f.end()); break;
            case 17: bii_normal.insert(bii_normal.end(), f.begin(), f.end()); break;
            case 18: aii_polar.insert(aii_polar.end(), f.begin(), f.end()); break;
            case 19: bii_polar.insert(bii_polar.end(), f.begin(), f.end()); break;
            case 20: aii14.insert(aii14.end(), f.begin(), f.end()); break;
            case 21: bii14.insert(bii14.end(), f.begin(), f.end()); break;
            case 23: ngbr14_flat += trim(line); break;
            case 24: {
                auto groups = checked_split_groups(f, 2);
                top_->ngbr14long.insert(top_->ngbr14long.end(), groups.begin(), groups.end());
                break;
            }
            case 25: ngbr23_flat += trim(line); break;
            case 26: {
                auto groups = checked_split_groups(f, 2);
                top_->ngbr23long.insert(top_->ngbr23long.end(), groups.begin(), groups.end());
                break;
            }
            case 29: top_->molecules.insert(top_->molecules.end(), f.begin(), f.end()); break;
            case 32:
                if (line.find("Exclusion") != std::string::npos && f.size() >= 2) {
                    if (parse_double(f[0]) > 30.0) throw parse_error("Sphere sizes exceeding 30A are currently not supported");
                    top_->exclusion = f[0];
                    top_->radii = f[1];
                }
                if ((line.find("Solute centre") != std::string::npos || line.find("Solute center") != std::string::npos) && f.size() >= 3) {
                    top_->solucenter = std::vector<std::string>(f.begin(), f.begin() + 3);
                }
                if ((line.find("Solvent centre") != std::string::npos || line.find("Solvent center") != std::string::npos) && f.size() >= 3) {
                    top_->solvcenter = std::vector<std::string>(f.begin(), f.begin() + 3);
                }
                break;
            case 33:
                for (char ch : line) {
                    if (!std::isspace(static_cast<unsigned char>(ch))) top_->excluded.push_back(ch == 'F' ? "0" : "1");
                }
                break;
            default: break;
        }
    }

    top_->coords = checked_split_groups(coord_flat, 3);
    top_->bonds = checked_split_groups(bond_flat, 3);
    top_->angles = checked_split_groups(angle_flat, 4);
    top_->torsions = checked_split_groups(torsion_flat, 5);
    top_->impropers = checked_split_groups(improper_flat, 5);
    for (size_t i = 0; i < ngbr14_flat.size(); i += line_width) top_->ngbr14.push_back(ngbr14_flat.substr(i, line_width));
    for (size_t i = 0; i < ngbr23_flat.size(); i += line_width) top_->ngbr23.push_back(ngbr23_flat.substr(i, line_width));

    std::map<std::string, int> charge_codes;
    int next_charge_code = 0;
    for (size_t i = 0; i < charges_tmp.size(); i++) {
        if (!charge_codes.count(charges_tmp[i])) charge_codes[charges_tmp[i]] = ++next_charge_code;
        top_->charges.push_back({static_cast<int>(i) + 1, charge_codes[charges_tmp[i]]});
    }
    for (const auto& item : charge_codes) {
        top_->ccharges[item.second] = item.first;
    }

    for (size_t i = 0; i < masses.size(); i++) {
        top_->catypes[static_cast<int>(i) + 1] = {
            masses[i],
            i < aii_normal.size() ? aii_normal[i] : "0",
            i < bii_normal.size() ? bii_normal[i] : "0",
            i < aii_polar.size() ? aii_polar[i] : "0",
            i < bii_polar.size() ? bii_polar[i] : "0",
            i < aii14.size() ? aii14[i] : "0",
            i < bii14.size() ? bii14[i] : "0",
        };
    }
    if (top_->solucenter.empty()) top_->solucenter = std::vector<std::string>(3, "0");
    if (top_->solvcenter.empty()) top_->solvcenter = std::vector<std::string>(3, "0");
}

void InpParser::ensure_fep() {
    if (fep_) return;
    ensure_md_input();
    ensure_topology();

    fep_ = std::make_unique<FepData>();
    if (fep_file_.empty()) return;

    std::ifstream first(fep_file_.c_str());
    if (!first) throw parse_error("Could not open FEP file " + fep_file_);
    std::string raw;
    while (std::getline(first, raw)) {
        std::vector<std::string> f = split_ws(strip_comment(raw));
        if (f.size() >= 2 && f[0] == "states") fep_->states = parse_int(f[1]);
    }

    fep_->q_atypes.assign(fep_->states, {});
    fep_->q_charges.assign(fep_->states, {});
    fep_->q_softpairs.assign(fep_->states, {});
    fep_->q_exclpairs.assign(fep_->states, {});
    fep_->q_softcores.assign(fep_->states, {});
    fep_->q_bonds.assign(fep_->states, {});
    fep_->q_angles.assign(fep_->states, {});
    fep_->q_torsions.assign(fep_->states, {});
    fep_->q_impropers.assign(fep_->states, {});
    fep_->q_shakes.assign(fep_->states, {});

    std::ifstream in(fep_file_.c_str());
    int block = 0;
    int atype_index = 0;
    while (std::getline(in, raw)) {
        std::string line = strip_comment(raw);
        if (line.empty()) continue;
        if (line.find("[atoms]") != std::string::npos) { block = 1; continue; }
        if (line.find("[FEP]") != std::string::npos) { block = 2; continue; }
        if (line.find("[change_charges]") != std::string::npos) { block = 3; continue; }
        if (line.find("[atom_types]") != std::string::npos) { block = 4; atype_index = 0; continue; }
        if (line.find("[change_atoms]") != std::string::npos) { block = 5; continue; }
        if (line.find("[soft_pairs]") != std::string::npos) { block = 6; continue; }
        if (line.find("[excluded_pairs]") != std::string::npos) { block = 7; continue; }
        if (line.find("[el_scale]") != std::string::npos) { block = 8; continue; }
        if (line.find("[softcore]") != std::string::npos) { block = 9; continue; }
        if (line.find("[bond_types]") != std::string::npos) { block = 12; fep_->q_cbonds.push_back({"0", "0.0", "0.0"}); continue; }
        if (line.find("[change_bonds]") != std::string::npos) { block = 13; continue; }
        if (line.find("[angle_types]") != std::string::npos) { block = 14; fep_->q_cangles.push_back({"0", "0.0", "0.0"}); continue; }
        if (line.find("[change_types]") != std::string::npos) { block = 15; continue; }
        if (line.find("[torsion_types]") != std::string::npos) { block = 16; fep_->q_ctorsions.push_back({"0", "0.0", "0.0", "0.0"}); continue; }
        if (line.find("[change_torsions]") != std::string::npos) { block = 17; continue; }
        if (line.find("[improper_types]") != std::string::npos) { block = 18; fep_->q_cimpropers.push_back({"0", "0.0", "0.0"}); continue; }
        if (line.find("[change_impropers]") != std::string::npos) { block = 19; continue; }
        if (line.find("[angle_couplings]") != std::string::npos) { block = 20; continue; }
        if (line.find("[torsion_couplings]") != std::string::npos) { block = 21; continue; }
        if (line.find("[improper_couplings]") != std::string::npos) { block = 22; continue; }
        if (line.find("[shake_constraints]") != std::string::npos) { block = 23; continue; }
        if (line.find("[off-diagonals]") != std::string::npos) { block = 24; continue; }

        std::vector<std::string> f = split_ws(line);
        if (f.empty()) continue;
        switch (block) {
            case 1:
                if (f.size() >= 2) fep_->q_atoms.push_back(f[1]);
                break;
            case 3:
                for (int s = 0; s < fep_->states && static_cast<size_t>(s + 1) < f.size(); s++) fep_->q_charges[s].push_back(f[s + 1]);
                break;
            case 4:
                atype_index++;
                fep_->q_catypes.push_back(f);
                fep_->atypemap[f[0]] = atype_index;
                break;
            case 5:
                for (int s = 0; s < fep_->states && static_cast<size_t>(s + 1) < f.size(); s++) {
                    auto it = fep_->atypemap.find(f[s + 1]);
                    fep_->q_atypes[s].push_back(it == fep_->atypemap.end() ? "0" : std::to_string(it->second));
                }
                break;
            case 6:
                for (int s = 0; s < fep_->states && static_cast<size_t>(s + 1) < f.size(); s++) fep_->q_softpairs[s].push_back(f[s + 1]);
                break;
            case 7:
                for (int s = 0; s < fep_->states && static_cast<size_t>(s + 1) < f.size(); s++) fep_->q_exclpairs[s].push_back(f[s + 1]);
                break;
            case 8: fep_->q_elscales.push_back(f); break;
            case 9:
                for (int s = 0; s < fep_->states && static_cast<size_t>(s + 1) < f.size(); s++) fep_->q_softcores[s].push_back(f[s + 1]);
                break;
            case 12: fep_->q_cbonds.push_back(f.size() > 1 ? std::vector<std::string>(f.begin() + 1, f.end()) : f); break;
            case 13:
                for (int s = 0; s < fep_->states && static_cast<size_t>(s + 1) < f.size(); s++) fep_->q_bonds[s].push_back(f[s + 1]);
                break;
            case 14: fep_->q_cangles.push_back(f.size() > 1 ? std::vector<std::string>(f.begin() + 1, f.end()) : f); break;
            case 15:
                for (int s = 0; s < fep_->states && static_cast<size_t>(s + 1) < f.size(); s++) fep_->q_angles[s].push_back(f[s + 1]);
                break;
            case 16: fep_->q_ctorsions.push_back(f.size() > 1 ? std::vector<std::string>(f.begin() + 1, f.end()) : f); break;
            case 17:
                for (int s = 0; s < fep_->states && static_cast<size_t>(s + 1) < f.size(); s++) fep_->q_torsions[s].push_back(f[s + 1]);
                break;
            case 18: fep_->q_cimpropers.push_back(f.size() > 1 ? std::vector<std::string>(f.begin() + 1, f.end()) : f); break;
            case 19:
                for (int s = 0; s < fep_->states && static_cast<size_t>(s + 1) < f.size(); s++) fep_->q_impropers[s].push_back(f[s + 1]);
                break;
            case 20: fep_->q_angcouples.push_back(f); break;
            case 21: fep_->q_torcouples.push_back(f); break;
            case 22: fep_->q_imprcouples.push_back(f); break;
            case 23:
                for (int s = 0; s < fep_->states && static_cast<size_t>(s + 1) < f.size(); s++) fep_->q_shakes[s].push_back(f[s + 1]);
                break;
            case 24: fep_->q_offdiags.push_back(f); break;
            default: break;
        }
    }

    if (fep_->states <= 0 || fep_->q_atoms.empty()) return;

    if (fep_->q_catypes.empty()) {
        for (const auto& item : top_->catypes) {
            const std::vector<std::string>& t = item.second;
            fep_->q_catypes.push_back({
                std::to_string(item.first),
                t.size() > 1 ? t[1] : "0",
                t.size() > 2 ? t[2] : "0",
                "0",
                "0",
                t.size() > 5 ? t[5] : "0",
                t.size() > 6 ? t[6] : "0",
                t.size() > 0 ? t[0] : "0",
            });
        }
    }

    for (int state = 0; state < fep_->states; state++) {
        if (fep_->q_atypes[state].empty()) {
            for (const std::string& q_atom : fep_->q_atoms) {
                int atom = parse_int(q_atom) - 1;
                if (atom < 0 || static_cast<size_t>(atom) >= top_->atypes.size()) throw parse_error("FEP atom index outside topology atom types.");
                fep_->q_atypes[state].push_back(std::to_string(top_->atypes[atom].second));
            }
        }

        if (fep_->q_charges[state].empty()) {
            for (const std::string& q_atom : fep_->q_atoms) {
                int atom = parse_int(q_atom) - 1;
                if (atom < 0 || static_cast<size_t>(atom) >= top_->charges.size()) throw parse_error("FEP atom index outside topology charges.");
                int charge_code = top_->charges[atom].second;
                auto charge = top_->ccharges.find(charge_code);
                if (charge == top_->ccharges.end()) throw parse_error("Topology atom charge code missing.");
                fep_->q_charges[state].push_back(charge->second);
            }
        }
    }
}

void InpParser::ensure_run_start_vectors() {
    if (run_start_ready_) return;
    ensure_md_input();
    ensure_topology();

    result.fresh_start = restart_file_.empty();
    if (!restart_file_.empty()) {
        read_restart_vectors(restart_file_, run_coords_, run_velocities_, result.restart_theta_corr);
    } else {
        run_coords_ = top_->coords;
        result.restart_theta_corr.clear();

        const auto& mdv = input_->keyed.count("md") ? input_->keyed["md"] : std::map<std::string, std::string>();
        int seed = parse_int(value_or(mdv, "random-seed", value_or(mdv, "random_seed", "1")));
        double temperature = parse_double(value_or(mdv, "initial-temperature", value_or(mdv, "initial_temperature", value_or(mdv, "temperature", "0"))));
        if (seed <= 0) throw parse_error("Fresh-start Q input requires a positive random_seed.");

        double kT = Boltz * temperature;
        run_velocities_.clear();
        run_velocities_.reserve(top_->atypes.size());
        for (size_t i = 0; i < top_->atypes.size(); i++) {
            int type_code = top_->atypes[i].second;
            auto catype = top_->catypes.find(type_code);
            if (catype == top_->catypes.end() || catype->second.empty()) throw parse_error("Topology atom type references missing mass.");
            double mass = parse_double(catype->second[0]);
            double sd = std::sqrt(kT / mass);
            run_velocities_.push_back({std::to_string(q_gauss(seed, sd)), std::to_string(q_gauss(seed, sd)), std::to_string(q_gauss(seed, sd))});
        }
    }

    if (run_coords_.size() != top_->coords.size() || run_velocities_.size() != top_->coords.size()) {
        throw parse_error("Input coordinate/velocity atom count does not match topology.");
    }
    run_start_ready_ = true;
}

void InpParser::parse_md() {
    ensure_md_input();
    const auto& files = input_->keyed["files"];
    const auto& mdv = input_->keyed.count("md") ? input_->keyed["md"] : std::map<std::string, std::string>();
    const auto& cut = input_->keyed.count("cut-offs") ? input_->keyed["cut-offs"] : std::map<std::string, std::string>();
    const auto& sphere = input_->keyed.count("sphere") ? input_->keyed["sphere"] : std::map<std::string, std::string>();
    const auto& solvent = input_->keyed.count("solvent") ? input_->keyed["solvent"] : std::map<std::string, std::string>();
    const auto& intervals = input_->keyed.count("intervals") ? input_->keyed["intervals"] : std::map<std::string, std::string>();

    md_t& md = result.md;
    md.steps = parse_int(value_or(mdv, "steps", "0"));
    md.stepsize = parse_double(value_or(mdv, "stepsize", "0"));
    md.temperature = parse_double(value_or(mdv, "temperature", "0"));
    md.thermostat = value_or(mdv, "thermostat", "berendsen");
    md.bath_coupling = parse_double(value_or(mdv, "bath-coupling", value_or(mdv, "bath_coupling", "1")));
    md.random_seed = parse_int(value_or(mdv, "random-seed", value_or(mdv, "random_seed", "1")));
    md.initial_temperature = parse_double(value_or(mdv, "initial-temperature", value_or(mdv, "initial_temperature", value_or(mdv, "temperature", "0"))));
    md.shake_solvent = is_on_value(mdv, "shake-solvent", bool_value(mdv, "shake_solvent", "off"));
    md.shake_solute = is_on_value(mdv, "shake-solute", bool_value(mdv, "shake_solute", "off"));
    md.shake_hydrogens = is_on_value(mdv, "shake-hydrogens", bool_value(mdv, "shake_hydrogens", "off"));
    md.lrf = is_on_value(mdv, "lrf", "off");
    md.separate_scaling = is_on_value(mdv, "separate-scaling", bool_value(mdv, "separate_scaling", "off"));
    md.charge_groups = true;
    md.solute_solute = parse_double(value_or(cut, "solute-solute", value_or(cut, "solute_solute", "10")));
    md.solvent_solvent = parse_double(value_or(cut, "solvent-solvent", value_or(cut, "solvent_solvent", "10")));
    md.solute_solvent = parse_double(value_or(cut, "solute-solvent", value_or(cut, "solute_solvent", "10")));
    md.q_atom = parse_double(value_or(cut, "q-atom", value_or(cut, "q_atom", "99")));
    md.shell_radius = parse_double(value_or(sphere, "shell-radius", value_or(sphere, "shell_radius", "0")));
    md.shell_force = parse_double(value_or(sphere, "shell-force", value_or(sphere, "shell_force", "10.0")));
    md.radial_force = parse_double(value_or(solvent, "radial-force", value_or(solvent, "radial_force", "60.0")));
    md.polarisation = true;
    md.polarisation_force = parse_double(value_or(solvent, "polarisation-force", value_or(solvent, "polarisation_force", "20.0")));
    md.non_bond = parse_int(value_or(intervals, "non-bond", value_or(intervals, "non_bond", "25")));
    md.output = parse_int(value_or(intervals, "output", "5"));
    md.energy = parse_int(value_or(intervals, "energy", "0"));
    md.trajectory = parse_int(value_or(intervals, "trajectory", "100"));

    result.native_output.final_file = resolve_output_path_local(value_or(files, "final", ""), input_dir_);
    result.native_output.trajectory_file = resolve_output_path_local(value_or(files, "trajectory", ""), input_dir_);
    result.native_output.energy_file = resolve_output_path_local(value_or(files, "energy", ""), input_dir_);
    result.native_output.trajectory_atoms = trajectory_atoms_setting(input_->lists);
    result.native_output.topology_file = topology_file_;
}

void InpParser::parse_lambdas() {
    ensure_md_input();
    result.lambdas.clear();
    auto lambdas_it = input_->lists.find("lambdas");
    if (lambdas_it != input_->lists.end()) {
        for (const auto& row : lambdas_it->second) {
            for (const std::string& value : row) result.lambdas.push_back(parse_double(value));
        }
    }

    auto read_rows = [&](const std::string& name) -> std::vector<std::vector<std::string>> {
        auto it = input_->lists.find(name);
        return it == input_->lists.end() ? std::vector<std::vector<std::string>>() : it->second;
    };

    auto seq = read_rows("sequence-restraints");
    result.restrseqs.resize(seq.size());
    for (size_t i = 0; i < seq.size(); i++) {
        result.restrseqs[i] = {row_int(seq[i], 0), row_int(seq[i], 1), row_double(seq[i], 2), row_int(seq[i], 3) == 1, row_int(seq[i], 4)};
    }

    auto pos = read_rows("atom-restraints");
    result.restrspos.resize(pos.size());
    for (size_t i = 0; i < pos.size(); i++) {
        result.restrspos[i] = {row_int(pos[i], 0), row_int(pos[i], 1), {row_double(pos[i], 2), row_double(pos[i], 3), row_double(pos[i], 4)}, {row_double(pos[i], 5), row_double(pos[i], 6), row_double(pos[i], 7)}};
    }

    auto dist = read_rows("distance-restraints");
    result.restrdists.resize(dist.size());
    for (size_t i = 0; i < dist.size(); i++) {
        restrdis_t r = {};
        r.ai = row_int(dist[i], 0);
        r.aj = row_int(dist[i], 1);
        r.d1 = row_double(dist[i], 2);
        r.d2 = row_double(dist[i], 3);
        r.k = row_double(dist[i], 4);
        r.ipsi = row_int(dist[i], 5);
        result.restrdists[i] = r;
    }

    auto angle = read_rows("angle-restraints");
    result.restrangs.resize(angle.size());
    for (size_t i = 0; i < angle.size(); i++) {
        result.restrangs[i] = {row_int(angle[i], 0), row_int(angle[i], 1), row_int(angle[i], 2), row_int(angle[i], 3), row_double(angle[i], 4), row_double(angle[i], 5)};
    }

    auto wall = read_rows("wall-restraints");
    result.restrwalls.resize(wall.size());
    for (size_t i = 0; i < wall.size(); i++) {
        result.restrwalls[i] = {row_int(wall[i], 0), row_int(wall[i], 1), row_double(wall[i], 2), row_double(wall[i], 3), row_double(wall[i], 5), row_double(wall[i], 4), row_int(wall[i], 6) == 1};
    }
}

void InpParser::parse_topo() {
    ensure_topology();
    result.topo.solvent_type = parse_int(top_->solvtype);
    result.topo.exclusion_radius = parse_double(top_->exclusion);
    result.topo.solvent_radius = parse_double(top_->radii);
    result.topo.solute_center = row_coord(top_->solucenter);
    result.topo.solvent_center = row_coord(top_->solvcenter);
    result.topo.el14_scale = parse_double(top_->scaling14);
    result.topo.coulomb_constant = parse_double(top_->coulomb);
    result.topo.vdw_rule = parse_int(top_->vdw_rule);
    if (result.topo.vdw_rule == 0) result.topo.vdw_rule = VDW_GEOMETRIC;
}

void InpParser::parse_coords_init() {
    ensure_topology();
    result.n_atoms_solute = parse_int(top_->natoms_solute);
    result.coords_init.resize(top_->coords.size());
    for (size_t i = 0; i < top_->coords.size(); i++) result.coords_init[i] = row_coord(top_->coords[i]);
}

void InpParser::parse_coords() {
    ensure_run_start_vectors();
    result.coords.resize(run_coords_.size());
    for (size_t i = 0; i < run_coords_.size(); i++) result.coords[i] = row_coord(run_coords_[i]);
}

void InpParser::parse_velocities() {
    ensure_run_start_vectors();
    result.velocities.resize(run_velocities_.size());
    for (size_t i = 0; i < run_velocities_.size(); i++) result.velocities[i] = row_vel(run_velocities_[i]);
}

void InpParser::parse_bonds() {
    ensure_topology();
    result.n_bonds_solute = parse_int(top_->nbonds_solute);
    result.bonds.resize(top_->bonds.size());
    for (size_t i = 0; i < top_->bonds.size(); i++) result.bonds[i] = {row_int(top_->bonds[i], 0), row_int(top_->bonds[i], 1), row_int(top_->bonds[i], 2)};
}

void InpParser::parse_cbonds() {
    ensure_topology();
    result.cbonds.clear();
    for (const auto& item : top_->cbonds) result.cbonds.push_back({item.first, row_double(item.second, 0), row_double(item.second, 1)});
}

void InpParser::parse_angles() {
    ensure_topology();
    result.n_angles_solute = parse_int(top_->nangles_solute);
    result.angles.resize(top_->angles.size());
    for (size_t i = 0; i < top_->angles.size(); i++) result.angles[i] = {row_int(top_->angles[i], 0), row_int(top_->angles[i], 1), row_int(top_->angles[i], 2), row_int(top_->angles[i], 3)};
}

void InpParser::parse_cangles() {
    ensure_topology();
    result.cangles.clear();
    for (const auto& item : top_->cangles) result.cangles.push_back({item.first, row_double(item.second, 0), row_double(item.second, 1)});
}

void InpParser::parse_torsions() {
    ensure_topology();
    result.n_torsions_solute = parse_int(top_->ntorsions_solute);
    result.torsions.resize(top_->torsions.size());
    for (size_t i = 0; i < top_->torsions.size(); i++) result.torsions[i] = {row_int(top_->torsions[i], 0), row_int(top_->torsions[i], 1), row_int(top_->torsions[i], 2), row_int(top_->torsions[i], 3), row_int(top_->torsions[i], 4)};
}

void InpParser::parse_ctorsions() {
    ensure_topology();
    result.ctorsions.clear();
    for (const auto& item : top_->ctorsions) {
        double paths = row_double(item.second, 3);
        result.ctorsions.push_back({item.first, row_double(item.second, 0), row_double(item.second, 1), row_double(item.second, 2), paths == 0.0 ? 0.0 : 1.0 / paths});
    }
}

void InpParser::parse_impropers() {
    ensure_topology();
    result.n_impropers_solute = parse_int(top_->nimpropers_solute);
    result.impropers.resize(top_->impropers.size());
    for (size_t i = 0; i < top_->impropers.size(); i++) result.impropers[i] = {row_int(top_->impropers[i], 0), row_int(top_->impropers[i], 1), row_int(top_->impropers[i], 2), row_int(top_->impropers[i], 3), row_int(top_->impropers[i], 4)};
}

void InpParser::parse_cimpropers() {
    ensure_topology();
    result.cimpropers.clear();
    for (const auto& item : top_->cimpropers) result.cimpropers.push_back({item.first, row_double(item.second, 0), row_double(item.second, 1)});
}

void InpParser::parse_restrspos() {}
void InpParser::parse_restrangs() {}
void InpParser::parse_restrdists() {}
void InpParser::parse_restrseqs() {}
void InpParser::parse_restrwalls() {}

void InpParser::parse_charges() {
    ensure_topology();
    result.charges.resize(top_->charges.size());
    for (size_t i = 0; i < top_->charges.size(); i++) result.charges[i] = {top_->charges[i].first, top_->charges[i].second};
}

void InpParser::parse_ccharges() {
    ensure_topology();
    result.ccharges.clear();
    for (const auto& item : top_->ccharges) result.ccharges.push_back({item.first, parse_double(item.second)});
}

void InpParser::parse_atypes() {
    ensure_topology();
    result.atypes.resize(top_->atypes.size());
    for (size_t i = 0; i < top_->atypes.size(); i++) result.atypes[i] = {top_->atypes[i].first, top_->atypes[i].second};
}

void InpParser::parse_catypes() {
    ensure_topology();
    result.catypes.clear();
    for (const auto& item : top_->catypes) {
        const auto& row = item.second;
        result.catypes.push_back({item.first, row_double(row, 0), row_double(row, 1), row_double(row, 2), row_double(row, 5), row_double(row, 6)});
    }
}

void InpParser::parse_heavy() {}

void InpParser::parse_excluded() {
    ensure_topology();
    result.excluded.assign(result.coords_init.size(), false);
    for (size_t i = 0; i < result.excluded.size() && i < top_->excluded.size(); i++) {
        result.excluded[i] = top_->excluded[i] == "1";
    }
}

void InpParser::parse_ngbrs14() {
    ensure_topology();
    add_short_neighbors(top_->ngbr14, result.ngbrs14);
}

void InpParser::parse_ngbrs14_long() {
    ensure_topology();
    add_long_neighbors(top_->ngbr14long, result.ngbrs14_long);
}

void InpParser::parse_ngbrs23() {
    ensure_topology();
    add_short_neighbors(top_->ngbr23, result.ngbrs23);
}

void InpParser::parse_ngbrs23_long() {
    ensure_topology();
    add_long_neighbors(top_->ngbr23long, result.ngbrs23_long);
}

void InpParser::parse_molecules() {
    ensure_topology();
    result.molecules.resize(top_->molecules.size());
    for (size_t i = 0; i < top_->molecules.size(); i++) result.molecules[i] = parse_int(top_->molecules[i]);
}

void InpParser::parse_charge_groups() {
    ensure_topology();
    charge_group_config_t config;
    config.n_cgrps_solute = parse_int(top_->solute_cgps);
    config.n_cgrps_solvent = parse_int(top_->solvent_cgps);
    config.iuse_switch_atom = 0;
    config.charge_groups.resize(top_->charge_group_headers.size());
    for (size_t i = 0; i < top_->charge_group_headers.size(); i++) {
        config.charge_groups[i].iswitch = row_int(top_->charge_group_headers[i], 1);
        config.charge_groups[i].atoms.resize(top_->charge_group_atoms[i].size());
        for (size_t j = 0; j < top_->charge_group_atoms[i].size(); j++) config.charge_groups[i].atoms[j] = parse_int(top_->charge_group_atoms[i][j]);
    }
    result.charge_groups.clear();
    result.charge_groups.push_back(config);
}

void InpParser::parse_q_atoms() {
    ensure_fep();
    result.q_atoms.resize(fep_->q_atoms.size());
    for (size_t i = 0; i < fep_->q_atoms.size(); i++) result.q_atoms[i] = parse_int(fep_->q_atoms[i]) - 1;
}

void InpParser::parse_q_angcouples() {
    ensure_fep();
    result.q_angcouples.resize(fep_->q_angcouples.size());
    for (size_t i = 0; i < fep_->q_angcouples.size(); i++) result.q_angcouples[i] = {row_int(fep_->q_angcouples[i], 0), row_int(fep_->q_angcouples[i], 1)};
}

void InpParser::parse_q_atypes() {
    ensure_fep();
    int n_lambdas = static_cast<int>(result.lambdas.size());
    result.q_atypes.resize(result.q_atoms.size() * n_lambdas);
    for (int state = 0; state < n_lambdas && state < static_cast<int>(fep_->q_atypes.size()); state++) {
        for (size_t i = 0; i < result.q_atoms.size() && i < fep_->q_atypes[state].size(); i++) {
            result.q_atypes[i + state * result.q_atoms.size()].code = parse_int(fep_->q_atypes[state][i]);
        }
    }
}

void InpParser::parse_q_catypes() {
    ensure_fep();
    result.q_catypes.resize(fep_->q_catypes.size());
    for (size_t i = 0; i < fep_->q_catypes.size(); i++) {
        const auto& row = fep_->q_catypes[i];
        result.q_catypes[i] = {static_cast<int>(i) + 1, row_double(row, 7), row_double(row, 1), row_double(row, 2), row_double(row, 5), row_double(row, 6)};
    }
}

void InpParser::parse_q_charges() {
    ensure_fep();
    int n_lambdas = static_cast<int>(result.lambdas.size());
    result.q_charges.resize(result.q_atoms.size() * n_lambdas);
    for (int state = 0; state < n_lambdas && state < static_cast<int>(fep_->q_charges.size()); state++) {
        for (size_t i = 0; i < result.q_atoms.size() && i < fep_->q_charges[state].size(); i++) {
            result.q_charges[i + state * result.q_atoms.size()].charge = parse_double(fep_->q_charges[state][i]);
        }
    }
}

void InpParser::parse_q_angles() {
    ensure_fep();
    int n_lambdas = static_cast<int>(result.lambdas.size());
    size_t total = 0;
    for (const auto& state : fep_->q_angles) total += state.size();
    int n_qangles = n_lambdas > 0 ? static_cast<int>(total) / n_lambdas : 0;
    result.q_angles.resize(n_qangles * n_lambdas);
    for (int state = 0; state < n_lambdas && state < static_cast<int>(fep_->q_angles.size()); state++) {
        for (int i = 0; i < n_qangles && static_cast<size_t>(i) < fep_->q_angles[state].size(); i++) {
            auto row = split_state_row(fep_->q_angles[state][i]);
            result.q_angles[i + state * n_qangles] = {row_int(row, 0), row_int(row, 1), row_int(row, 2), row_int(row, 3)};
        }
    }
}

void InpParser::parse_q_cangles() {
    ensure_fep();
    result.q_cangles.resize(fep_->q_cangles.size());
    for (size_t i = 0; i < fep_->q_cangles.size(); i++) result.q_cangles[i] = {0, row_double(fep_->q_cangles[i], 0), row_double(fep_->q_cangles[i], 1)};
}

void InpParser::parse_q_bonds() {
    ensure_fep();
    int n_lambdas = static_cast<int>(result.lambdas.size());
    size_t total = 0;
    for (const auto& state : fep_->q_bonds) total += state.size();
    int n_qbonds = n_lambdas > 0 ? static_cast<int>(total) / n_lambdas : 0;
    result.q_bonds.resize(n_qbonds * n_lambdas);
    for (int state = 0; state < n_lambdas && state < static_cast<int>(fep_->q_bonds.size()); state++) {
        for (int i = 0; i < n_qbonds && static_cast<size_t>(i) < fep_->q_bonds[state].size(); i++) {
            auto row = split_state_row(fep_->q_bonds[state][i]);
            result.q_bonds[i + state * n_qbonds] = {row_int(row, 0), row_int(row, 1), row_int(row, 2)};
        }
    }
}

void InpParser::parse_q_cbonds() {
    ensure_fep();
    result.q_cbonds.resize(fep_->q_cbonds.size());
    for (size_t i = 0; i < fep_->q_cbonds.size(); i++) result.q_cbonds[i] = {0, row_double(fep_->q_cbonds[i], 0), row_double(fep_->q_cbonds[i], 1)};
}

void InpParser::parse_q_impropers() {
    ensure_fep();
    int n_lambdas = static_cast<int>(result.lambdas.size());
    size_t total = 0;
    for (const auto& state : fep_->q_impropers) total += state.size();
    int n_qimpropers = n_lambdas > 0 ? static_cast<int>(total) / n_lambdas : 0;
    result.q_impropers.resize(n_qimpropers * n_lambdas);
    for (int state = 0; state < n_lambdas && state < static_cast<int>(fep_->q_impropers.size()); state++) {
        for (int i = 0; i < n_qimpropers && static_cast<size_t>(i) < fep_->q_impropers[state].size(); i++) {
            auto row = split_state_row(fep_->q_impropers[state][i]);
            result.q_impropers[i + state * n_qimpropers] = {row_int(row, 0), row_int(row, 1), row_int(row, 2), row_int(row, 3), row_int(row, 4)};
        }
    }
}

void InpParser::parse_q_cimpropers() {
    ensure_fep();
    result.q_cimpropers.resize(fep_->q_cimpropers.size());
    for (size_t i = 0; i < fep_->q_cimpropers.size(); i++) result.q_cimpropers[i] = {row_double(fep_->q_cimpropers[i], 0), row_double(fep_->q_cimpropers[i], 1)};
}

void InpParser::parse_q_torsions() {
    ensure_fep();
    int n_lambdas = static_cast<int>(result.lambdas.size());
    size_t total = 0;
    for (const auto& state : fep_->q_torsions) total += state.size();
    int n_qtorsions = n_lambdas > 0 ? static_cast<int>(total) / n_lambdas : 0;
    result.q_torsions.resize(n_qtorsions * n_lambdas);
    for (int state = 0; state < n_lambdas && state < static_cast<int>(fep_->q_torsions.size()); state++) {
        for (int i = 0; i < n_qtorsions && static_cast<size_t>(i) < fep_->q_torsions[state].size(); i++) {
            auto row = split_state_row(fep_->q_torsions[state][i]);
            result.q_torsions[i + state * n_qtorsions] = {row_int(row, 0), row_int(row, 1), row_int(row, 2), row_int(row, 3), row_int(row, 4)};
        }
    }
}

void InpParser::parse_q_ctorsions() {
    ensure_fep();
    result.q_ctorsions.resize(fep_->q_ctorsions.size());
    for (size_t i = 0; i < fep_->q_ctorsions.size(); i++) result.q_ctorsions[i] = {0, row_double(fep_->q_ctorsions[i], 0), row_double(fep_->q_ctorsions[i], 1), row_double(fep_->q_ctorsions[i], 2), 0.0};
}

void InpParser::parse_q_offdiags() {
    ensure_fep();
    result.q_offdiags.resize(fep_->q_offdiags.size());
    for (size_t i = 0; i < fep_->q_offdiags.size(); i++) {
        result.q_offdiags[i] = {row_int(fep_->q_offdiags[i], 0), row_int(fep_->q_offdiags[i], 1), row_int(fep_->q_offdiags[i], 2), row_int(fep_->q_offdiags[i], 3), row_double(fep_->q_offdiags[i], 4), row_double(fep_->q_offdiags[i], 5)};
    }
}

void InpParser::parse_q_imprcouples() {
    ensure_fep();
    result.q_imprcouples.resize(fep_->q_imprcouples.size());
    for (size_t i = 0; i < fep_->q_imprcouples.size(); i++) result.q_imprcouples[i] = {row_int(fep_->q_imprcouples[i], 0), row_int(fep_->q_imprcouples[i], 1)};
}

void InpParser::parse_q_softpairs() {
    ensure_fep();
    result.q_softpairs.clear();
    for (const auto& state : fep_->q_softpairs) {
        for (const std::string& value : state) {
            auto row = split_state_row(value);
            result.q_softpairs.push_back({row_int(row, 0), row_int(row, 1)});
        }
    }
}

void InpParser::parse_q_torcouples() {
    ensure_fep();
    result.q_torcouples.resize(fep_->q_torcouples.size());
    for (size_t i = 0; i < fep_->q_torcouples.size(); i++) result.q_torcouples[i] = {row_int(fep_->q_torcouples[i], 0), row_int(fep_->q_torcouples[i], 1)};
}

void InpParser::parse_q_elscales() {
    ensure_fep();
    result.q_elscales.resize(fep_->q_elscales.size());
    for (size_t i = 0; i < fep_->q_elscales.size(); i++) result.q_elscales[i] = {row_int(fep_->q_elscales[i], 0), row_int(fep_->q_elscales[i], 1), row_double(fep_->q_elscales[i], 2)};
}

void InpParser::parse_q_exclpairs() {
    ensure_fep();
    int n_lambdas = static_cast<int>(result.lambdas.size());
    size_t total = 0;
    for (const auto& state : fep_->q_exclpairs) total += state.size();
    int n_qexclpairs = n_lambdas > 0 ? static_cast<int>(total) / n_lambdas : 0;
    result.q_exclpairs.resize(n_qexclpairs * n_lambdas);
    for (int state = 0; state < n_lambdas && state < static_cast<int>(fep_->q_exclpairs.size()); state++) {
        for (int i = 0; i < n_qexclpairs && static_cast<size_t>(i) < fep_->q_exclpairs[state].size(); i++) {
            auto row = split_state_row(fep_->q_exclpairs[state][i]);
            result.q_exclpairs[i + state * n_qexclpairs] = {row_int(row, 0), row_int(row, 1), row_int(row, 2)};
        }
    }
}

void InpParser::parse_q_shakes() {
    ensure_fep();
    int n_lambdas = static_cast<int>(result.lambdas.size());
    size_t total = 0;
    for (const auto& state : fep_->q_shakes) total += state.size();
    int n_qshakes = n_lambdas > 0 ? static_cast<int>(total) / n_lambdas : 0;
    result.q_shakes.resize(n_qshakes * n_lambdas);
    for (int state = 0; state < n_lambdas && state < static_cast<int>(fep_->q_shakes.size()); state++) {
        for (int i = 0; i < n_qshakes && static_cast<size_t>(i) < fep_->q_shakes[state].size(); i++) {
            auto row = split_state_row(fep_->q_shakes[state][i]);
            result.q_shakes[i + state * n_qshakes] = {row_int(row, 0), row_int(row, 1), row_double(row, 2)};
        }
    }
}

void InpParser::parse_q_softcores() {
    ensure_fep();
    int n_lambdas = static_cast<int>(result.lambdas.size());
    size_t total = 0;
    for (const auto& state : fep_->q_softcores) total += state.size();
    int n_qsoftcores = n_lambdas > 0 ? static_cast<int>(total) / n_lambdas : 0;
    result.q_softcores.resize(n_qsoftcores * n_lambdas);
    for (int state = 0; state < n_lambdas && state < static_cast<int>(fep_->q_softcores.size()); state++) {
        for (int i = 0; i < n_qsoftcores && static_cast<size_t>(i) < fep_->q_softcores[state].size(); i++) {
            result.q_softcores[i + state * n_qsoftcores].s = parse_double(fep_->q_softcores[state][i]);
        }
    }
}
