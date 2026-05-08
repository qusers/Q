#include "q_input.h"
#include "q_output.h"

#include "constants.h"
#include "context.h"

#include <algorithm>
#include <cctype>
#include <cmath>
#include <cstdint>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <fstream>
#include <map>
#include <memory>
#include <sstream>
#include <string>
#include <unistd.h>
#include <vector>

q_input_mode_t q_input_mode = Q_INPUT_CSV;
bool q_native_fresh_start = false;

static void fatal(const std::string &message) {
    printf(">>> FATAL: %s\n", message.c_str());
    exit(EXIT_FAILURE);
}

static bool file_exists(const std::string &path) {
    return access(path.c_str(), F_OK) == 0;
}

static std::string trim(const std::string &value) {
    size_t first = value.find_first_not_of(" \t\r\n");
    if (first == std::string::npos) return "";
    size_t last = value.find_last_not_of(" \t\r\n");
    return value.substr(first, last - first + 1);
}

static std::string strip_comment(const std::string &line) {
    size_t bang = line.find('!');
    if (bang == std::string::npos) return trim(line);
    return trim(line.substr(0, bang));
}

static std::string lower_normalized(std::string value) {
    std::transform(value.begin(), value.end(), value.begin(), [](unsigned char c) {
        if (c == '_') return '-';
        return (char) std::tolower(c);
    });
    return value;
}

static std::vector<std::string> split_ws(const std::string &line) {
    std::vector<std::string> fields;
    std::istringstream iss(line);
    std::string field;
    while (iss >> field) fields.push_back(field);
    return fields;
}

static std::vector<std::vector<std::string> > split_groups(const std::vector<std::string> &flat, size_t width) {
    std::vector<std::vector<std::string> > groups;
    for (size_t i = 0; i + width <= flat.size(); i += width) {
        groups.push_back(std::vector<std::string>(flat.begin() + i, flat.begin() + i + width));
    }
    return groups;
}

static std::string dirname_of(const std::string &path) {
    size_t slash = path.find_last_of('/');
    if (slash == std::string::npos) return ".";
    if (slash == 0) return "/";
    return path.substr(0, slash);
}

static bool is_absolute(const std::string &path) {
    return !path.empty() && path[0] == '/';
}

static std::string join_path(const std::string &left, const std::string &right) {
    if (left.empty() || left == ".") return right;
    if (!right.empty() && right[0] == '/') return right;
    if (left[left.size() - 1] == '/') return left + right;
    return left + "/" + right;
}

static std::string resolve_existing_or_cwd(const std::string &value, const std::string &input_dir) {
    if (value.empty()) return "";
    if (is_absolute(value)) return value;
    std::string input_relative = join_path(input_dir, value);
    if (file_exists(input_relative)) return input_relative;
    char cwd[1024];
    if (getcwd(cwd, sizeof(cwd)) == NULL) fatal("Could not determine current working directory.");
    return join_path(cwd, value);
}

static std::string resolve_output_path(const std::string &value, const std::string &input_dir) {
    if (value.empty()) return "";
    if (is_absolute(value)) return value;
    return join_path(input_dir, value);
}

static std::string bool_value(const std::map<std::string, std::string> &values, const std::string &key, const std::string &fallback) {
    std::map<std::string, std::string>::const_iterator it = values.find(key);
    if (it == values.end()) return fallback;
    std::string value = lower_normalized(it->second);
    if (value == "true" || value == "yes" || value == "1") return "on";
    if (value == "false" || value == "no" || value == "0") return "off";
    return value;
}

static std::string value_or(const std::map<std::string, std::string> &values, const std::string &key, const std::string &fallback) {
    std::map<std::string, std::string>::const_iterator it = values.find(key);
    if (it == values.end()) return fallback;
    return it->second;
}

struct MDInput {
    std::map<std::string, std::map<std::string, std::string> > keyed;
    std::map<std::string, std::vector<std::vector<std::string> > > lists;
};

static MDInput parse_md_input(const std::string &path) {
    std::ifstream in(path.c_str());
    if (!in) fatal("Could not open Q input file " + path);

    MDInput input;
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
            input.lists[section].push_back(fields);
        }
        else if (fields.size() >= 2) {
            std::string key = lower_normalized(fields[0]);
            std::string value = fields[1];
            for (size_t i = 2; i < fields.size(); i++) value += " " + fields[i];
            input.keyed[section][key] = value;
        }
    }
    return input;
}

struct ChargeGroup {
    std::vector<std::string> header;
    std::vector<std::string> atoms;
};

struct TopData {
    std::string natoms_solute = "0";
    std::string nbonds_solute = "0";
    std::string nangles_solute = "0";
    std::string ntorsions_solute = "0";
    std::string nimpropers_solute = "0";
    std::vector<std::vector<std::string> > coords, bonds, angles, torsions, impropers, ngbr14long, ngbr23long;
    std::vector<std::pair<int, int> > atypes;
    std::map<int, std::vector<std::string> > cbonds, cangles, ctorsions, cimpropers, catypes;
    std::vector<std::pair<int, int> > charges;
    std::map<int, std::string> ccharges;
    std::vector<std::string> ngbr14, ngbr23, molecules, excluded;
    std::string solute_cgps = "0", solvent_cgps = "0";
    std::vector<ChargeGroup> charge_groups;
    std::string solvtype = "0", exclusion = "0", radii = "0", scaling14 = "1", coulomb = "332.0716", vdw_rule = "1";
    std::vector<std::string> solucenter, solvcenter;
};

static TopData parse_topology(const std::string &path) {
    std::ifstream in(path.c_str());
    if (!in) fatal("Could not open topology file " + path);

    TopData data;
    int block = 0;
    int atype_count = 0;
    int charge_group_switch = 0;
    int charge_group_atoms_expected = 0;
    ChargeGroup current_group;
    std::vector<std::string> coord_flat, bond_flat, angle_flat, torsion_flat, improper_flat;
    std::vector<std::string> charges_tmp, masses, aii_normal, bii_normal, aii_polar, bii_polar, aii14, bii14;
    std::string ngbr14_flat, ngbr23_flat;

    std::string raw;
    while (std::getline(in, raw)) {
        std::string line = trim(raw);
        if (line.empty()) continue;

        if (line.find("Total no. of atoms") != std::string::npos) {
            std::vector<std::string> f = split_ws(line);
            if (f.size() > 1) data.natoms_solute = f[1];
            block = 1;
            continue;
        }
        if (line.find("No. of integer atom codes") != std::string::npos) { block = 2; continue; }
        if (line.find("No. of bonds") != std::string::npos) {
            std::vector<std::string> f = split_ws(line);
            if (f.size() > 1) data.nbonds_solute = f[1];
            block = 3;
            continue;
        }
        if (line.find("No. of bond codes") != std::string::npos) { block = 4; continue; }
        if (line.find("No. of angles") != std::string::npos) {
            std::vector<std::string> f = split_ws(line);
            if (f.size() > 1) data.nangles_solute = f[1];
            block = 5;
            continue;
        }
        if (line.find("No. of angle codes") != std::string::npos) { block = 6; continue; }
        if (line.find("No. of torsions") != std::string::npos) {
            std::vector<std::string> f = split_ws(line);
            if (f.size() > 1) data.ntorsions_solute = f[1];
            block = 7;
            continue;
        }
        if (line.find("No. of torsion codes") != std::string::npos) { block = 8; continue; }
        if (line.find("No. of impropers") != std::string::npos) {
            std::vector<std::string> f = split_ws(line);
            if (f.size() > 1) data.nimpropers_solute = f[1];
            block = 9;
            continue;
        }
        if (line.find("No. of improper codes") != std::string::npos) { block = 10; continue; }
        if (line.find("No. of atomic charges") != std::string::npos) { block = 11; continue; }
        if (line.find("No. of charge groups") != std::string::npos) {
            std::vector<std::string> f = split_ws(line);
            int total = f.size() > 0 ? atoi(f[0].c_str()) : 0;
            int solute = f.size() > 1 ? atoi(f[1].c_str()) : 0;
            data.solute_cgps = std::to_string(solute);
            data.solvent_cgps = std::to_string(total - solute);
            block = 12;
            charge_group_switch = 1;
            continue;
        }
        if (line.find("vdW combination rule") != std::string::npos) {
            std::vector<std::string> f = split_ws(line);
            if (!f.empty()) data.vdw_rule = f[0];
            block = 13;
            continue;
        }
        if (line.find("Electrostatic 1-4 scaling factor") != std::string::npos) { block = 14; }
        if (line.find("Masses") != std::string::npos) { block = 15; continue; }
        if (line.find("sqrt (Aii) normal") != std::string::npos) { block = 16; continue; }
        if (line.find("sqrt (Bii) normal") != std::string::npos) { block = 17; continue; }
        if (line.find("sqrt (Aii) polar") != std::string::npos) { block = 18; continue; }
        if (line.find("sqrt (Bii) polar") != std::string::npos) { block = 19; continue; }
        if (line.find("sqrt (Aii) 1-4") != std::string::npos) { block = 20; continue; }
        if (line.find("sqrt (Bii) 1-4") != std::string::npos) { block = 21; continue; }
        if (line.find("R* normal:") != std::string::npos) { block = 16; continue; }
        if (line.find("epsilon normal:") != std::string::npos) { block = 17; continue; }
        if (line.find("R* polar:") != std::string::npos) { block = 18; continue; }
        if (line.find("epsilon polar:") != std::string::npos) { block = 19; continue; }
        if (line.find("R* 1-4:") != std::string::npos) { block = 20; continue; }
        if (line.find("epsilon 1-4:") != std::string::npos) { block = 21; continue; }
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
            if (!f.empty()) data.solvtype = f[0];
            if (data.solvtype != "0") fatal("solvent type other than SPC (0) is not supported");
            block = 32;
            continue;
        }
        if (line.find("No. of excluded atoms") != std::string::npos) { block = 33; continue; }

        std::vector<std::string> f = split_ws(line);
        switch (block) {
            case 1: coord_flat.insert(coord_flat.end(), f.begin(), f.end()); break;
            case 2:
                for (size_t i = 0; i < f.size(); i++) data.atypes.push_back(std::make_pair(++atype_count, atoi(f[i].c_str())));
                break;
            case 3: bond_flat.insert(bond_flat.end(), f.begin(), f.end()); break;
            case 4:
                if (f.size() >= 3) data.cbonds[atoi(f[0].c_str())] = std::vector<std::string>(f.begin() + 1, f.begin() + 3);
                break;
            case 5: angle_flat.insert(angle_flat.end(), f.begin(), f.end()); break;
            case 6:
                if (f.size() >= 3) data.cangles[atoi(f[0].c_str())] = std::vector<std::string>(f.begin() + 1, f.begin() + 3);
                break;
            case 7: torsion_flat.insert(torsion_flat.end(), f.begin(), f.end()); break;
            case 8:
                if (f.size() >= 5) data.ctorsions[atoi(f[0].c_str())] = std::vector<std::string>(f.begin() + 1, f.begin() + 5);
                break;
            case 9: improper_flat.insert(improper_flat.end(), f.begin(), f.end()); break;
            case 10:
                if (f.size() >= 3) data.cimpropers[atoi(f[0].c_str())] = std::vector<std::string>(f.begin() + 1, f.begin() + 3);
                break;
            case 11: charges_tmp.insert(charges_tmp.end(), f.begin(), f.end()); break;
            case 12:
                if (charge_group_switch == 1 && f.size() >= 2) {
                    current_group = ChargeGroup();
                    current_group.header = f;
                    charge_group_atoms_expected = atoi(f[0].c_str());
                    charge_group_switch = 2;
                }
                else if (charge_group_switch == 2) {
                    current_group.atoms.insert(current_group.atoms.end(), f.begin(), f.end());
                    if ((int) current_group.atoms.size() >= charge_group_atoms_expected) {
                        data.charge_groups.push_back(current_group);
                        charge_group_switch = 1;
                    }
                }
                break;
            case 14:
                if (f.size() >= 2) {
                    data.scaling14 = f[0];
                    data.coulomb = f[1];
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
                std::vector<std::vector<std::string> > groups = split_groups(f, 2);
                data.ngbr14long.insert(data.ngbr14long.end(), groups.begin(), groups.end());
                break;
            }
            case 25: ngbr23_flat += trim(line); break;
            case 26: {
                std::vector<std::vector<std::string> > groups = split_groups(f, 2);
                data.ngbr23long.insert(data.ngbr23long.end(), groups.begin(), groups.end());
                break;
            }
            case 29: data.molecules.insert(data.molecules.end(), f.begin(), f.end()); break;
            case 32:
                if (line.find("Exclusion") != std::string::npos && f.size() >= 2) {
                    if (atof(f[0].c_str()) > 30.0) fatal("Sphere sizes exceeding 30A are currently not supported");
                    data.exclusion = f[0];
                    data.radii = f[1];
                }
                if ((line.find("Solute centre") != std::string::npos || line.find("Solute center") != std::string::npos) && f.size() >= 3) {
                    data.solucenter = std::vector<std::string>(f.begin(), f.begin() + 3);
                }
                if ((line.find("Solvent centre") != std::string::npos || line.find("Solvent center") != std::string::npos) && f.size() >= 3) {
                    data.solvcenter = std::vector<std::string>(f.begin(), f.begin() + 3);
                }
                break;
            case 33:
                for (size_t i = 0; i < line.size(); i++) {
                    if (!isspace((unsigned char) line[i])) data.excluded.push_back(line[i] == 'F' ? "0" : "1");
                }
                break;
            default: break;
        }
    }

    data.coords = split_groups(coord_flat, 3);
    data.bonds = split_groups(bond_flat, 3);
    data.angles = split_groups(angle_flat, 4);
    data.torsions = split_groups(torsion_flat, 5);
    data.impropers = split_groups(improper_flat, 5);
    for (size_t i = 0; i < ngbr14_flat.size(); i += 25) data.ngbr14.push_back(ngbr14_flat.substr(i, 25));
    for (size_t i = 0; i < ngbr23_flat.size(); i += 25) data.ngbr23.push_back(ngbr23_flat.substr(i, 25));

    std::map<std::string, int> charge_codes;
    int next_charge_code = 0;
    for (size_t i = 0; i < charges_tmp.size(); i++) {
        if (!charge_codes.count(charges_tmp[i])) charge_codes[charges_tmp[i]] = ++next_charge_code;
        data.charges.push_back(std::make_pair((int) i + 1, charge_codes[charges_tmp[i]]));
    }
    for (std::map<std::string, int>::iterator it = charge_codes.begin(); it != charge_codes.end(); ++it) {
        data.ccharges[it->second] = it->first;
    }

    size_t ntypes = masses.size();
    for (size_t i = 0; i < ntypes; i++) {
        std::vector<std::string> row;
        row.push_back(masses[i]);
        row.push_back(i < aii_normal.size() ? aii_normal[i] : "0");
        row.push_back(i < bii_normal.size() ? bii_normal[i] : "0");
        row.push_back(i < aii_polar.size() ? aii_polar[i] : "0");
        row.push_back(i < bii_polar.size() ? bii_polar[i] : "0");
        row.push_back(i < aii14.size() ? aii14[i] : "0");
        row.push_back(i < bii14.size() ? bii14[i] : "0");
        data.catypes[(int) i + 1] = row;
    }
    if (data.solucenter.empty()) data.solucenter = std::vector<std::string>(3, "0");
    if (data.solvcenter.empty()) data.solvcenter = std::vector<std::string>(3, "0");
    return data;
}

struct FepData {
    int states = 0;
    std::vector<std::string> q_atoms;
    std::vector<std::vector<std::string> > q_atypes, q_charges, q_softpairs, q_exclpairs, q_softcores, q_bonds, q_angles, q_torsions, q_impropers, q_shakes;
    std::vector<std::vector<std::string> > q_catypes, q_elscales, q_cbonds, q_cangles, q_ctorsions, q_cimpropers, q_angcouples, q_torcouples, q_imprcouples, q_offdiags;
    std::map<std::string, int> atypemap;
};

static void ensure_state_vectors(FepData &data) {
    data.q_atypes.assign(data.states, std::vector<std::string>());
    data.q_charges.assign(data.states, std::vector<std::string>());
    data.q_softpairs.assign(data.states, std::vector<std::string>());
    data.q_exclpairs.assign(data.states, std::vector<std::string>());
    data.q_softcores.assign(data.states, std::vector<std::string>());
    data.q_bonds.assign(data.states, std::vector<std::string>());
    data.q_angles.assign(data.states, std::vector<std::string>());
    data.q_torsions.assign(data.states, std::vector<std::string>());
    data.q_impropers.assign(data.states, std::vector<std::string>());
    data.q_shakes.assign(data.states, std::vector<std::string>());
}

static FepData parse_fep(const std::string &path) {
    FepData data;
    if (path.empty()) return data;

    std::ifstream first(path.c_str());
    if (!first) fatal("Could not open FEP file " + path);
    std::string raw;
    while (std::getline(first, raw)) {
        std::vector<std::string> f = split_ws(strip_comment(raw));
        if (f.size() >= 2 && f[0] == "states") data.states = atoi(f[1].c_str());
    }
    ensure_state_vectors(data);

    std::ifstream in(path.c_str());
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
        if (line.find("[bond_types]") != std::string::npos) { block = 12; data.q_cbonds.push_back(std::vector<std::string>{"0", "0.0", "0.0"}); continue; }
        if (line.find("[change_bonds]") != std::string::npos) { block = 13; continue; }
        if (line.find("[angle_types]") != std::string::npos) { block = 14; data.q_cangles.push_back(std::vector<std::string>{"0", "0.0", "0.0"}); continue; }
        if (line.find("[change_types]") != std::string::npos) { block = 15; continue; }
        if (line.find("[torsion_types]") != std::string::npos) { block = 16; data.q_ctorsions.push_back(std::vector<std::string>{"0", "0.0", "0.0", "0.0"}); continue; }
        if (line.find("[change_torsions]") != std::string::npos) { block = 17; continue; }
        if (line.find("[improper_types]") != std::string::npos) { block = 18; data.q_cimpropers.push_back(std::vector<std::string>{"0", "0.0", "0.0"}); continue; }
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
                if (f.size() >= 2) data.q_atoms.push_back(f[1]);
                break;
            case 3:
                for (int s = 0; s < data.states && (size_t) s + 1 < f.size(); s++) data.q_charges[s].push_back(f[s + 1]);
                break;
            case 4:
                atype_index++;
                data.q_catypes.push_back(f);
                data.atypemap[f[0]] = atype_index;
                break;
            case 5:
                for (int s = 0; s < data.states && (size_t) s + 1 < f.size(); s++) {
                    std::map<std::string, int>::iterator it = data.atypemap.find(f[s + 1]);
                    data.q_atypes[s].push_back(it == data.atypemap.end() ? "0" : std::to_string(it->second));
                }
                break;
            case 6:
                for (int s = 0; s < data.states && (size_t) s + 1 < f.size(); s++) data.q_softpairs[s].push_back(f[s + 1]);
                break;
            case 7:
                for (int s = 0; s < data.states && (size_t) s + 1 < f.size(); s++) data.q_exclpairs[s].push_back(f[s + 1]);
                break;
            case 8: data.q_elscales.push_back(f); break;
            case 9:
                for (int s = 0; s < data.states && (size_t) s + 1 < f.size(); s++) data.q_softcores[s].push_back(f[s + 1]);
                break;
            case 12: data.q_cbonds.push_back(f.size() > 1 ? std::vector<std::string>(f.begin() + 1, f.end()) : f); break;
            case 13:
                for (int s = 0; s < data.states && (size_t) s + 1 < f.size(); s++) data.q_bonds[s].push_back(f[s + 1]);
                break;
            case 14: data.q_cangles.push_back(f.size() > 1 ? std::vector<std::string>(f.begin() + 1, f.end()) : f); break;
            case 15:
                for (int s = 0; s < data.states && (size_t) s + 1 < f.size(); s++) data.q_angles[s].push_back(f[s + 1]);
                break;
            case 16: data.q_ctorsions.push_back(f.size() > 1 ? std::vector<std::string>(f.begin() + 1, f.end()) : f); break;
            case 17:
                for (int s = 0; s < data.states && (size_t) s + 1 < f.size(); s++) data.q_torsions[s].push_back(f[s + 1]);
                break;
            case 18: data.q_cimpropers.push_back(f.size() > 1 ? std::vector<std::string>(f.begin() + 1, f.end()) : f); break;
            case 19:
                for (int s = 0; s < data.states && (size_t) s + 1 < f.size(); s++) data.q_impropers[s].push_back(f[s + 1]);
                break;
            case 20: data.q_angcouples.push_back(f); break;
            case 21: data.q_torcouples.push_back(f); break;
            case 22: data.q_imprcouples.push_back(f); break;
            case 23:
                for (int s = 0; s < data.states && (size_t) s + 1 < f.size(); s++) data.q_shakes[s].push_back(f[s + 1]);
                break;
            case 24: data.q_offdiags.push_back(f); break;
            default: break;
        }
    }
    return data;
}

static void complete_fep_defaults(FepData &fep, const TopData &top) {
    if (fep.states <= 0 || fep.q_atoms.empty()) return;

    if (fep.q_catypes.empty()) {
        for (std::map<int, std::vector<std::string> >::const_iterator it = top.catypes.begin(); it != top.catypes.end(); ++it) {
            const std::vector<std::string> &t = it->second;
            std::vector<std::string> row;
            row.push_back(std::to_string(it->first));
            row.push_back(t.size() > 1 ? t[1] : "0");
            row.push_back(t.size() > 2 ? t[2] : "0");
            row.push_back("0");
            row.push_back("0");
            row.push_back(t.size() > 5 ? t[5] : "0");
            row.push_back(t.size() > 6 ? t[6] : "0");
            row.push_back(t.size() > 0 ? t[0] : "0");
            fep.q_catypes.push_back(row);
        }
    }

    for (int state = 0; state < fep.states; state++) {
        if ((size_t) state >= fep.q_atypes.size() || (size_t) state >= fep.q_charges.size()) {
            fatal("Invalid FEP state table allocation.");
        }

        if (fep.q_atypes[state].empty()) {
            for (size_t i = 0; i < fep.q_atoms.size(); i++) {
                int atom = atoi(fep.q_atoms[i].c_str()) - 1;
                if (atom < 0 || (size_t) atom >= top.atypes.size()) fatal("FEP atom index outside topology atom types.");
                fep.q_atypes[state].push_back(std::to_string(top.atypes[atom].second));
            }
        }

        if (fep.q_charges[state].empty()) {
            for (size_t i = 0; i < fep.q_atoms.size(); i++) {
                int atom = atoi(fep.q_atoms[i].c_str()) - 1;
                if (atom < 0 || (size_t) atom >= top.charges.size()) fatal("FEP atom index outside topology charges.");
                int charge_code = top.charges[atom].second;
                std::map<int, std::string>::const_iterator charge = top.ccharges.find(charge_code);
                if (charge == top.ccharges.end()) fatal("Topology atom charge code missing.");
                fep.q_charges[state].push_back(charge->second);
            }
        }
    }
}

static int to_int(const std::string &value) {
    return atoi(value.c_str());
}

static double to_double(const std::string &value) {
    return atof(value.c_str());
}

static int row_int(const std::vector<std::string> &row, size_t index, int fallback = 0) {
    return index < row.size() ? to_int(row[index]) : fallback;
}

static double row_double(const std::vector<std::string> &row, size_t index, double fallback = 0.0) {
    return index < row.size() ? to_double(row[index]) : fallback;
}

static bool is_on_value(const std::map<std::string, std::string> &values, const std::string &key, const std::string &fallback) {
    return bool_value(values, key, fallback) == "on";
}

static void set_md_context(const MDInput &input) {
    Context &ctx = Context::instance();
    std::map<std::string, std::string> mdv = input.keyed.count("md") ? input.keyed.find("md")->second : std::map<std::string, std::string>();
    std::map<std::string, std::string> cut = input.keyed.count("cut-offs") ? input.keyed.find("cut-offs")->second : std::map<std::string, std::string>();
    std::map<std::string, std::string> sphere = input.keyed.count("sphere") ? input.keyed.find("sphere")->second : std::map<std::string, std::string>();
    std::map<std::string, std::string> solvent = input.keyed.count("solvent") ? input.keyed.find("solvent")->second : std::map<std::string, std::string>();
    std::map<std::string, std::string> intervals = input.keyed.count("intervals") ? input.keyed.find("intervals")->second : std::map<std::string, std::string>();

    md_t &md = ctx.md;
    md.steps = atoi(value_or(mdv, "steps", "0").c_str());
    md.stepsize = atof(value_or(mdv, "stepsize", "0").c_str());
    md.temperature = atof(value_or(mdv, "temperature", "0").c_str());
    std::string thermostat = value_or(mdv, "thermostat", "berendsen");
    strncpy(md.thermostat, thermostat.c_str(), sizeof(md.thermostat) - 1);
    md.thermostat[sizeof(md.thermostat) - 1] = '\0';
    md.bath_coupling = atof(value_or(mdv, "bath-coupling", value_or(mdv, "bath_coupling", "1")).c_str());
    md.random_seed = atoi(value_or(mdv, "random-seed", value_or(mdv, "random_seed", "1")).c_str());
    md.initial_temperature = atof(value_or(mdv, "initial-temperature", value_or(mdv, "initial_temperature", value_or(mdv, "temperature", "0"))).c_str());
    md.shake_solvent = is_on_value(mdv, "shake-solvent", bool_value(mdv, "shake_solvent", "off"));
    md.shake_solute = is_on_value(mdv, "shake-solute", bool_value(mdv, "shake_solute", "off"));
    md.shake_hydrogens = is_on_value(mdv, "shake-hydrogens", bool_value(mdv, "shake_hydrogens", "off"));
    md.lrf = is_on_value(mdv, "lrf", "off");
    md.separate_scaling = is_on_value(mdv, "separate-scaling", bool_value(mdv, "separate_scaling", "off"));
    md.charge_groups = true;
    md.solute_solute = atof(value_or(cut, "solute-solute", value_or(cut, "solute_solute", "10")).c_str());
    md.solvent_solvent = atof(value_or(cut, "solvent-solvent", value_or(cut, "solvent_solvent", "10")).c_str());
    md.solute_solvent = atof(value_or(cut, "solute-solvent", value_or(cut, "solute_solvent", "10")).c_str());
    md.q_atom = atof(value_or(cut, "q-atom", value_or(cut, "q_atom", "99")).c_str());
    md.shell_radius = atof(value_or(sphere, "shell-radius", value_or(sphere, "shell_radius", "0")).c_str());
    md.shell_force = atof(value_or(sphere, "shell-force", value_or(sphere, "shell_force", "10.0")).c_str());
    md.radial_force = atof(value_or(solvent, "radial-force", value_or(solvent, "radial_force", "60.0")).c_str());
    md.polarisation = true;
    md.polarisation_force = atof(value_or(solvent, "polarisation-force", value_or(solvent, "polarisation_force", "20.0")).c_str());
    md.non_bond = atoi(value_or(intervals, "non-bond", value_or(intervals, "non_bond", "25")).c_str());
    md.output = atoi(value_or(intervals, "output", "5").c_str());
    md.energy = atoi(value_or(intervals, "energy", "1").c_str());
    md.trajectory = atoi(value_or(intervals, "trajectory", "100").c_str());

    std::vector<std::string> lambdas_flat;
    std::map<std::string, std::vector<std::vector<std::string> > >::const_iterator lambdas_it = input.lists.find("lambdas");
    if (lambdas_it != input.lists.end()) {
        for (size_t i = 0; i < lambdas_it->second.size(); i++) {
            lambdas_flat.insert(lambdas_flat.end(), lambdas_it->second[i].begin(), lambdas_it->second[i].end());
        }
    }
    ctx.n_lambdas = (int) lambdas_flat.size();
    ctx.lambdas = std::make_unique<HostDeviceBuffer<double>>(ctx.n_lambdas, true, ctx.run_gpu);
    for (int i = 0; i < ctx.n_lambdas; i++) ctx.lambdas->cpu_data_p[i] = to_double(lambdas_flat[i]);

    std::vector<std::vector<std::string> > seq = input.lists.count("sequence-restraints") ? input.lists.find("sequence-restraints")->second : std::vector<std::vector<std::string> >();
    ctx.n_restrseqs = (int) seq.size();
    ctx.restrseqs = std::make_unique<HostDeviceBuffer<restrseq_t>>(ctx.n_restrseqs, true, ctx.run_gpu);
    for (int i = 0; i < ctx.n_restrseqs; i++) {
        restrseq_t r = {};
        r.ai = row_int(seq[i], 0);
        r.aj = row_int(seq[i], 1);
        r.k = row_double(seq[i], 2);
        r.ih = row_int(seq[i], 3) == 1;
        r.to_center = row_int(seq[i], 4);
        ctx.restrseqs->cpu_data_p[i] = r;
    }

    std::vector<std::vector<std::string> > pos = input.lists.count("atom-restraints") ? input.lists.find("atom-restraints")->second : std::vector<std::vector<std::string> >();
    ctx.n_restrspos = (int) pos.size();
    ctx.restrspos = std::make_unique<HostDeviceBuffer<restrpos_t>>(ctx.n_restrspos, true, ctx.run_gpu);
    for (int i = 0; i < ctx.n_restrspos; i++) {
        restrpos_t r = {};
        r.a = row_int(pos[i], 0);
        r.ipsi = row_int(pos[i], 1);
        r.x = {row_double(pos[i], 2), row_double(pos[i], 3), row_double(pos[i], 4)};
        r.k = {row_double(pos[i], 5), row_double(pos[i], 6), row_double(pos[i], 7)};
        ctx.restrspos->cpu_data_p[i] = r;
    }

    std::vector<std::vector<std::string> > dist = input.lists.count("distance-restraints") ? input.lists.find("distance-restraints")->second : std::vector<std::vector<std::string> >();
    ctx.n_restrdists = (int) dist.size();
    ctx.restrdists = std::make_unique<HostDeviceBuffer<restrdis_t>>(ctx.n_restrdists, true, ctx.run_gpu);
    for (int i = 0; i < ctx.n_restrdists; i++) {
        restrdis_t r = {};
        r.ai = row_int(dist[i], 0);
        r.aj = row_int(dist[i], 1);
        r.d1 = row_double(dist[i], 2);
        r.d2 = row_double(dist[i], 3);
        r.k = row_double(dist[i], 4);
        r.ipsi = row_int(dist[i], 5);
        ctx.restrdists->cpu_data_p[i] = r;
    }

    std::vector<std::vector<std::string> > angle = input.lists.count("angle-restraints") ? input.lists.find("angle-restraints")->second : std::vector<std::vector<std::string> >();
    ctx.n_restrangs = (int) angle.size();
    ctx.restrangs = std::make_unique<HostDeviceBuffer<restrang_t>>(ctx.n_restrangs, true, ctx.run_gpu);
    for (int i = 0; i < ctx.n_restrangs; i++) {
        restrang_t r = {};
        r.ai = row_int(angle[i], 0);
        r.aj = row_int(angle[i], 1);
        r.ak = row_int(angle[i], 2);
        r.ipsi = row_int(angle[i], 3);
        r.ang = row_double(angle[i], 4);
        r.k = row_double(angle[i], 5);
        ctx.restrangs->cpu_data_p[i] = r;
    }

    std::vector<std::vector<std::string> > wall = input.lists.count("wall-restraints") ? input.lists.find("wall-restraints")->second : std::vector<std::vector<std::string> >();
    ctx.n_restrwalls = (int) wall.size();
    ctx.restrwalls = std::make_unique<HostDeviceBuffer<restrwall_t>>(ctx.n_restrwalls, true, ctx.run_gpu);
    for (int i = 0; i < ctx.n_restrwalls; i++) {
        restrwall_t r = {};
        r.ai = row_int(wall[i], 0);
        r.aj = row_int(wall[i], 1);
        r.d = row_double(wall[i], 2);
        r.k = row_double(wall[i], 3);
        r.dMorse = row_double(wall[i], 4);
        r.aMorse = row_double(wall[i], 5);
        r.ih = row_int(wall[i], 6) == 1;
        ctx.restrwalls->cpu_data_p[i] = r;
    }
}

static void set_lj_from_short_lines(const std::vector<std::string> &lines, int value) {
    Context &ctx = Context::instance();
    int *matrix = ctx.LJ_matrix->cpu_data_p;
    int n_lines = (int) lines.size();
    for (int line_i = 0; line_i < n_lines; line_i++) {
        for (int i = 0; i < line_width && (size_t) i < lines[line_i].size(); i++) {
            if (lines[line_i][i] == '1') {
                int j = (line_i + i + 1) % n_lines;
                matrix[line_i * ctx.n_atoms_solute + j] = value;
                matrix[j * ctx.n_atoms_solute + line_i] = value;
            }
        }
    }
}

static void set_lj_from_long_rows(const std::vector<std::vector<std::string> > &rows, int value) {
    Context &ctx = Context::instance();
    int *matrix = ctx.LJ_matrix->cpu_data_p;
    for (size_t i = 0; i < rows.size(); i++) {
        int ai = row_int(rows[i], 0) - 1;
        int aj = row_int(rows[i], 1) - 1;
        if (ai >= 0 && aj >= 0 && ai < ctx.n_atoms_solute && aj < ctx.n_atoms_solute) {
            matrix[ai * ctx.n_atoms_solute + aj] = value;
            matrix[aj * ctx.n_atoms_solute + ai] = value;
        }
    }
}

static coord_t row_coord(const std::vector<std::string> &row) {
    return {row_double(row, 0), row_double(row, 1), row_double(row, 2)};
}

static void set_topology_context(const TopData &top, const std::vector<std::vector<std::string> > &input_coords) {
    Context &ctx = Context::instance();

    ctx.n_atoms = (int) top.coords.size();
    ctx.n_atoms_solute = to_int(top.natoms_solute);
    ctx.coords_init = std::make_unique<HostDeviceBuffer<coord_t>>(ctx.n_atoms, true, ctx.run_gpu);
    ctx.coords = std::make_unique<HostDeviceBuffer<coord_t>>(ctx.n_atoms, true, ctx.run_gpu);
    for (int i = 0; i < ctx.n_atoms; i++) {
        ctx.coords_init->cpu_data_p[i] = row_coord(top.coords[i]);
        ctx.coords->cpu_data_p[i] = row_coord(input_coords[i]);
    }

    ctx.topo.solvent_type = to_int(top.solvtype);
    ctx.topo.exclusion_radius = to_double(top.exclusion);
    ctx.topo.solvent_radius = to_double(top.radii);
    ctx.topo.solute_center = row_coord(top.solucenter);
    ctx.topo.solvent_center = row_coord(top.solvcenter);
    ctx.topo.el14_scale = to_double(top.scaling14);
    ctx.topo.coulomb_constant = to_double(top.coulomb);
    ctx.topo.vdw_rule = to_int(top.vdw_rule);
    if (ctx.topo.vdw_rule == 0) ctx.topo.vdw_rule = VDW_GEOMETRIC;

    ctx.n_atypes = (int) top.atypes.size();
    ctx.atypes = std::make_unique<HostDeviceBuffer<atype_t>>(ctx.n_atypes, true, ctx.run_gpu);
    for (int i = 0; i < ctx.n_atypes; i++) ctx.atypes->cpu_data_p[i] = {top.atypes[i].first, top.atypes[i].second};

    ctx.n_catypes = (int) top.catypes.size();
    ctx.catypes = std::make_unique<HostDeviceBuffer<catype_t>>(ctx.n_catypes, true, ctx.run_gpu);
    int idx = 0;
    for (std::map<int, std::vector<std::string> >::const_iterator it = top.catypes.begin(); it != top.catypes.end(); ++it, ++idx) {
        const std::vector<std::string> &r = it->second;
        ctx.catypes->cpu_data_p[idx] = {it->first, row_double(r, 0), row_double(r, 1), row_double(r, 2), row_double(r, 5), row_double(r, 6)};
    }
    if (ctx.topo.vdw_rule == VDW_ARITHMETIC) {
        for (int i = 0; i < ctx.n_catypes; i++) {
            ctx.catypes->cpu_data_p[i].bii_normal = sqrt(fabs(ctx.catypes->cpu_data_p[i].bii_normal));
            ctx.catypes->cpu_data_p[i].bii_1_4 = sqrt(fabs(ctx.catypes->cpu_data_p[i].bii_1_4));
        }
    }

    ctx.n_charges = (int) top.charges.size();
    ctx.charges = std::make_unique<HostDeviceBuffer<charge_t>>(ctx.n_charges, true, ctx.run_gpu);
    for (int i = 0; i < ctx.n_charges; i++) ctx.charges->cpu_data_p[i] = {top.charges[i].first, top.charges[i].second};

    ctx.n_ccharges = (int) top.ccharges.size();
    ctx.ccharges = std::make_unique<HostDeviceBuffer<ccharge_t>>(ctx.n_ccharges, true, ctx.run_gpu);
    idx = 0;
    for (std::map<int, std::string>::const_iterator it = top.ccharges.begin(); it != top.ccharges.end(); ++it, ++idx) {
        ctx.ccharges->cpu_data_p[idx] = {it->first, to_double(it->second)};
    }

    ctx.n_bonds = (int) top.bonds.size();
    ctx.n_bonds_solute = to_int(top.nbonds_solute);
    ctx.bonds = std::make_unique<HostDeviceBuffer<bond_t>>(ctx.n_bonds, true, ctx.run_gpu);
    for (int i = 0; i < ctx.n_bonds; i++) ctx.bonds->cpu_data_p[i] = {row_int(top.bonds[i], 0), row_int(top.bonds[i], 1), row_int(top.bonds[i], 2)};

    ctx.n_cbonds = (int) top.cbonds.size();
    ctx.cbonds = std::make_unique<HostDeviceBuffer<cbond_t>>(ctx.n_cbonds, true, ctx.run_gpu);
    idx = 0;
    for (std::map<int, std::vector<std::string> >::const_iterator it = top.cbonds.begin(); it != top.cbonds.end(); ++it, ++idx) {
        ctx.cbonds->cpu_data_p[idx] = {it->first, row_double(it->second, 0), row_double(it->second, 1)};
    }

    ctx.n_angles = (int) top.angles.size();
    ctx.n_angles_solute = to_int(top.nangles_solute);
    ctx.angles = std::make_unique<HostDeviceBuffer<angle_t>>(ctx.n_angles, true, ctx.run_gpu);
    for (int i = 0; i < ctx.n_angles; i++) ctx.angles->cpu_data_p[i] = {row_int(top.angles[i], 0), row_int(top.angles[i], 1), row_int(top.angles[i], 2), row_int(top.angles[i], 3)};

    ctx.n_cangles = (int) top.cangles.size();
    ctx.cangles = std::make_unique<HostDeviceBuffer<cangle_t>>(ctx.n_cangles, true, ctx.run_gpu);
    idx = 0;
    for (std::map<int, std::vector<std::string> >::const_iterator it = top.cangles.begin(); it != top.cangles.end(); ++it, ++idx) {
        ctx.cangles->cpu_data_p[idx] = {it->first, row_double(it->second, 0), row_double(it->second, 1)};
    }

    ctx.n_torsions = (int) top.torsions.size();
    ctx.n_torsions_solute = to_int(top.ntorsions_solute);
    ctx.torsions = std::make_unique<HostDeviceBuffer<torsion_t>>(ctx.n_torsions, true, ctx.run_gpu);
    for (int i = 0; i < ctx.n_torsions; i++) ctx.torsions->cpu_data_p[i] = {row_int(top.torsions[i], 0), row_int(top.torsions[i], 1), row_int(top.torsions[i], 2), row_int(top.torsions[i], 3), row_int(top.torsions[i], 4)};

    ctx.n_ctorsions = (int) top.ctorsions.size();
    ctx.ctorsions = std::make_unique<HostDeviceBuffer<ctorsion_t>>(ctx.n_ctorsions, true, ctx.run_gpu);
    idx = 0;
    for (std::map<int, std::vector<std::string> >::const_iterator it = top.ctorsions.begin(); it != top.ctorsions.end(); ++it, ++idx) {
        double paths = row_double(it->second, 3);
        ctx.ctorsions->cpu_data_p[idx] = {it->first, row_double(it->second, 0), row_double(it->second, 1), row_double(it->second, 2), paths == 0.0 ? 0.0 : 1.0 / paths};
    }

    ctx.n_impropers = (int) top.impropers.size();
    ctx.n_impropers_solute = to_int(top.nimpropers_solute);
    ctx.impropers = std::make_unique<HostDeviceBuffer<improper_t>>(ctx.n_impropers, true, ctx.run_gpu);
    for (int i = 0; i < ctx.n_impropers; i++) ctx.impropers->cpu_data_p[i] = {row_int(top.impropers[i], 0), row_int(top.impropers[i], 1), row_int(top.impropers[i], 2), row_int(top.impropers[i], 3), row_int(top.impropers[i], 4)};

    ctx.n_cimpropers = (int) top.cimpropers.size();
    ctx.cimpropers = std::make_unique<HostDeviceBuffer<cimproper_t>>(ctx.n_cimpropers, true, ctx.run_gpu);
    idx = 0;
    for (std::map<int, std::vector<std::string> >::const_iterator it = top.cimpropers.begin(); it != top.cimpropers.end(); ++it, ++idx) {
        ctx.cimpropers->cpu_data_p[idx] = {it->first, row_double(it->second, 0), row_double(it->second, 1)};
    }

    ctx.n_molecules = (int) top.molecules.size();
    ctx.molecules.resize(ctx.n_molecules);
    for (int i = 0; i < ctx.n_molecules; i++) ctx.molecules[i] = to_int(top.molecules[i]);

    ctx.excluded = std::make_unique<HostDeviceBuffer<bool>>(ctx.n_atoms, true, ctx.run_gpu);
    ctx.n_excluded = 0;
    for (int i = 0; i < ctx.n_atoms; i++) {
        bool excl = (size_t) i < top.excluded.size() && top.excluded[i] == "1";
        ctx.excluded->cpu_data_p[i] = excl;
        if (excl) ctx.n_excluded++;
    }

    ctx.LJ_matrix = std::make_unique<HostDeviceBuffer<int>>(ctx.n_atoms_solute * ctx.n_atoms_solute, true, ctx.run_gpu);
    set_lj_from_short_lines(top.ngbr14, 1);
    set_lj_from_short_lines(top.ngbr23, 3);
    set_lj_from_long_rows(top.ngbr14long, 1);
    set_lj_from_long_rows(top.ngbr23long, 3);

    charge_group_config_t config;
    config.n_cgrps_solute = to_int(top.solute_cgps);
    config.n_cgrps_solvent = to_int(top.solvent_cgps);
    config.iuse_switch_atom = 0;
    config.charge_groups.resize(top.charge_groups.size());
    for (size_t i = 0; i < top.charge_groups.size(); i++) {
        charge_group_t cg;
        cg.iswitch = row_int(top.charge_groups[i].header, 1);
        cg.atoms.resize(top.charge_groups[i].atoms.size());
        for (size_t j = 0; j < top.charge_groups[i].atoms.size(); j++) cg.atoms[j] = to_int(top.charge_groups[i].atoms[j]);
        config.charge_groups[i] = cg;
    }
    ctx.charge_group_config = config;
}

static void set_velocities_context(const std::vector<std::vector<std::string> > &input_velocities) {
    Context &ctx = Context::instance();
    ctx.velocities = std::make_unique<HostDeviceBuffer<vel_t>>(ctx.n_atoms, true, ctx.run_gpu);
    for (int i = 0; i < ctx.n_atoms; i++) {
        ctx.velocities->cpu_data_p[i] = {row_double(input_velocities[i], 0), row_double(input_velocities[i], 1), row_double(input_velocities[i], 2)};
    }
}

static void set_flat_q_bonds(const std::vector<std::vector<std::string> > &state_rows) {
    Context &ctx = Context::instance();
    size_t count = 0;
    for (size_t s = 0; s < state_rows.size(); s++) count += state_rows[s].size();
    ctx.n_qbonds = ctx.n_lambdas > 0 ? (int) count / ctx.n_lambdas : 0;
    ctx.q_bonds.resize(ctx.n_qbonds * ctx.n_lambdas);
    int k = 0;
    for (int s = 0; s < ctx.n_lambdas && (size_t) s < state_rows.size(); s++) {
        for (size_t i = 0; i < state_rows[s].size() && i < (size_t) ctx.n_qbonds; i++) {
            std::vector<std::string> row = split_ws(state_rows[s][i]);
            ctx.q_bonds[(int) i + s * ctx.n_qbonds] = {row_int(row, 0), row_int(row, 1), row_int(row, 2)};
            k++;
        }
    }
    (void) k;
}

static void set_fep_context(const FepData &fep) {
    Context &ctx = Context::instance();
    if (fep.states > 0 && fep.states != ctx.n_lambdas) {
        fatal("FEP state count must match [lambdas] count in native QGPU input.");
    }

    ctx.n_qatoms = (int) fep.q_atoms.size();
    ctx.q_atoms.resize(ctx.n_qatoms);
    for (int i = 0; i < ctx.n_qatoms; i++) ctx.q_atoms[i] = to_int(fep.q_atoms[i]) - 1;

    ctx.n_qangcouples = (int) fep.q_angcouples.size();
    ctx.q_angcouples.resize(ctx.n_qangcouples);
    for (int i = 0; i < ctx.n_qangcouples; i++) ctx.q_angcouples[i] = {row_int(fep.q_angcouples[i], 0), row_int(fep.q_angcouples[i], 1)};

    ctx.n_qcatypes = (int) fep.q_catypes.size();
    ctx.q_catypes.resize(ctx.n_qcatypes);
    for (int i = 0; i < ctx.n_qcatypes; i++) {
        const std::vector<std::string> &r = fep.q_catypes[i];
        ctx.q_catypes[i].code = i + 1;
        ctx.q_catypes[i].aii_normal = row_double(r, 1);
        ctx.q_catypes[i].bii_normal = row_double(r, 2);
        ctx.q_catypes[i].aii_1_4 = row_double(r, 5);
        ctx.q_catypes[i].bii_1_4 = row_double(r, 6);
        ctx.q_catypes[i].m = row_double(r, 7);
    }
    if (ctx.topo.vdw_rule == VDW_ARITHMETIC) {
        for (int i = 0; i < ctx.n_qcatypes; i++) {
            ctx.q_catypes[i].bii_normal = sqrt(fabs(ctx.q_catypes[i].bii_normal));
            ctx.q_catypes[i].bii_1_4 = sqrt(fabs(ctx.q_catypes[i].bii_1_4));
        }
    }

    ctx.q_atypes.resize(ctx.n_qatoms * ctx.n_lambdas);
    for (int s = 0; s < ctx.n_lambdas && (size_t) s < fep.q_atypes.size(); s++) {
        for (int i = 0; i < ctx.n_qatoms && (size_t) i < fep.q_atypes[s].size(); i++) {
            ctx.q_atypes[i + s * ctx.n_qatoms].code = to_int(fep.q_atypes[s][i]);
        }
    }

    ctx.q_charges.resize(ctx.n_qatoms * ctx.n_lambdas);
    for (int s = 0; s < ctx.n_lambdas && (size_t) s < fep.q_charges.size(); s++) {
        for (int i = 0; i < ctx.n_qatoms && (size_t) i < fep.q_charges[s].size(); i++) {
            ctx.q_charges[i + s * ctx.n_qatoms].charge = to_double(fep.q_charges[s][i]);
        }
    }

    ctx.n_qcangles = (int) fep.q_cangles.size();
    ctx.q_cangles.resize(ctx.n_qcangles);
    for (int i = 0; i < ctx.n_qcangles; i++) {
        ctx.q_cangles[i].kth = row_double(fep.q_cangles[i], 0);
        ctx.q_cangles[i].th0 = row_double(fep.q_cangles[i], 1);
    }

    ctx.n_qcbonds = (int) fep.q_cbonds.size();
    ctx.q_cbonds.resize(ctx.n_qcbonds);
    for (int i = 0; i < ctx.n_qcbonds; i++) {
        ctx.q_cbonds[i].kb = row_double(fep.q_cbonds[i], 0);
        ctx.q_cbonds[i].b0 = row_double(fep.q_cbonds[i], 1);
    }

    ctx.n_qcimpropers = (int) fep.q_cimpropers.size();
    ctx.q_cimpropers.resize(ctx.n_qcimpropers);
    for (int i = 0; i < ctx.n_qcimpropers; i++) ctx.q_cimpropers[i] = {row_double(fep.q_cimpropers[i], 0), row_double(fep.q_cimpropers[i], 1)};

    ctx.n_qctorsions = (int) fep.q_ctorsions.size();
    ctx.q_ctorsions.resize(ctx.n_qctorsions);
    for (int i = 0; i < ctx.n_qctorsions; i++) {
        ctx.q_ctorsions[i].k = row_double(fep.q_ctorsions[i], 0);
        ctx.q_ctorsions[i].n = row_double(fep.q_ctorsions[i], 1);
        ctx.q_ctorsions[i].d = row_double(fep.q_ctorsions[i], 2);
    }

    ctx.n_qoffdiags = (int) fep.q_offdiags.size();
    ctx.q_offdiags.resize(ctx.n_qoffdiags);
    for (int i = 0; i < ctx.n_qoffdiags; i++) {
        ctx.q_offdiags[i] = {row_int(fep.q_offdiags[i], 0), row_int(fep.q_offdiags[i], 1), row_int(fep.q_offdiags[i], 2), row_int(fep.q_offdiags[i], 3),
                             row_double(fep.q_offdiags[i], 4), row_double(fep.q_offdiags[i], 5)};
    }

    ctx.n_qimprcouples = (int) fep.q_imprcouples.size();
    ctx.q_imprcouples.resize(ctx.n_qimprcouples);
    for (int i = 0; i < ctx.n_qimprcouples; i++) ctx.q_imprcouples[i] = {row_int(fep.q_imprcouples[i], 0), row_int(fep.q_imprcouples[i], 1)};

    ctx.n_qsoftpairs = 0;
    for (size_t s = 0; s < fep.q_softpairs.size(); s++) ctx.n_qsoftpairs += (int) fep.q_softpairs[s].size();
    ctx.q_softpairs.resize(ctx.n_qsoftpairs);
    int out = 0;
    for (size_t s = 0; s < fep.q_softpairs.size(); s++) {
        for (size_t i = 0; i < fep.q_softpairs[s].size(); i++) {
            std::vector<std::string> row = split_ws(fep.q_softpairs[s][i]);
            ctx.q_softpairs[out++] = {row_int(row, 0), row_int(row, 1)};
        }
    }

    ctx.n_qtorcouples = (int) fep.q_torcouples.size();
    ctx.q_torcouples.resize(ctx.n_qtorcouples);
    for (int i = 0; i < ctx.n_qtorcouples; i++) ctx.q_torcouples[i] = {row_int(fep.q_torcouples[i], 0), row_int(fep.q_torcouples[i], 1)};

    set_flat_q_bonds(fep.q_bonds);

    size_t q_angle_count = 0;
    for (size_t s = 0; s < fep.q_angles.size(); s++) q_angle_count += fep.q_angles[s].size();
    ctx.n_qangles = ctx.n_lambdas > 0 ? (int) q_angle_count / ctx.n_lambdas : 0;
    ctx.q_angles.resize(ctx.n_qangles * ctx.n_lambdas);
    for (int s = 0; s < ctx.n_lambdas && (size_t) s < fep.q_angles.size(); s++) {
        for (int i = 0; i < ctx.n_qangles && (size_t) i < fep.q_angles[s].size(); i++) {
            std::vector<std::string> row = split_ws(fep.q_angles[s][i]);
            ctx.q_angles[i + s * ctx.n_qangles] = {row_int(row, 0), row_int(row, 1), row_int(row, 2), row_int(row, 3)};
        }
    }

    size_t q_improper_count = 0;
    for (size_t s = 0; s < fep.q_impropers.size(); s++) q_improper_count += fep.q_impropers[s].size();
    ctx.n_qimpropers = ctx.n_lambdas > 0 ? (int) q_improper_count / ctx.n_lambdas : 0;
    ctx.q_impropers.resize(ctx.n_qimpropers * ctx.n_lambdas);
    for (int s = 0; s < ctx.n_lambdas && (size_t) s < fep.q_impropers.size(); s++) {
        for (int i = 0; i < ctx.n_qimpropers && (size_t) i < fep.q_impropers[s].size(); i++) {
            std::vector<std::string> row = split_ws(fep.q_impropers[s][i]);
            ctx.q_impropers[i + s * ctx.n_qimpropers] = {row_int(row, 0), row_int(row, 1), row_int(row, 2), row_int(row, 3), row_int(row, 4)};
        }
    }

    size_t q_torsion_count = 0;
    for (size_t s = 0; s < fep.q_torsions.size(); s++) q_torsion_count += fep.q_torsions[s].size();
    ctx.n_qtorsions = ctx.n_lambdas > 0 ? (int) q_torsion_count / ctx.n_lambdas : 0;
    ctx.q_torsions.resize(ctx.n_qtorsions * ctx.n_lambdas);
    for (int s = 0; s < ctx.n_lambdas && (size_t) s < fep.q_torsions.size(); s++) {
        for (int i = 0; i < ctx.n_qtorsions && (size_t) i < fep.q_torsions[s].size(); i++) {
            std::vector<std::string> row = split_ws(fep.q_torsions[s][i]);
            ctx.q_torsions[i + s * ctx.n_qtorsions] = {row_int(row, 0), row_int(row, 1), row_int(row, 2), row_int(row, 3), row_int(row, 4)};
        }
    }

    size_t q_shake_count = 0;
    for (size_t s = 0; s < fep.q_shakes.size(); s++) q_shake_count += fep.q_shakes[s].size();
    ctx.n_qshakes = ctx.n_lambdas > 0 ? (int) q_shake_count / ctx.n_lambdas : 0;
    ctx.q_shakes.resize(ctx.n_qshakes * ctx.n_lambdas);
    for (int s = 0; s < ctx.n_lambdas && (size_t) s < fep.q_shakes.size(); s++) {
        for (int i = 0; i < ctx.n_qshakes && (size_t) i < fep.q_shakes[s].size(); i++) {
            std::vector<std::string> row = split_ws(fep.q_shakes[s][i]);
            ctx.q_shakes[i + s * ctx.n_qshakes] = {row_int(row, 0), row_int(row, 1), row_double(row, 2)};
        }
    }

    size_t q_softcore_count = 0;
    for (size_t s = 0; s < fep.q_softcores.size(); s++) q_softcore_count += fep.q_softcores[s].size();
    ctx.n_qsoftcores = ctx.n_lambdas > 0 ? (int) q_softcore_count / ctx.n_lambdas : 0;
    ctx.q_softcores.resize(ctx.n_qsoftcores * ctx.n_lambdas);
    for (int s = 0; s < ctx.n_lambdas && (size_t) s < fep.q_softcores.size(); s++) {
        for (int i = 0; i < ctx.n_qsoftcores && (size_t) i < fep.q_softcores[s].size(); i++) {
            ctx.q_softcores[i + s * ctx.n_qsoftcores].s = to_double(fep.q_softcores[s][i]);
        }
    }

    ctx.n_qelscales = ctx.n_lambdas > 0 ? (int) fep.q_elscales.size() / ctx.n_lambdas : 0;
    ctx.q_elscales = std::make_unique<HostDeviceBuffer<q_elscale_t>>(ctx.n_qelscales * ctx.n_lambdas, true, ctx.run_gpu);
    for (int s = 0; s < ctx.n_lambdas; s++) {
        for (int i = 0; i < ctx.n_qelscales; i++) {
            const std::vector<std::string> &row = fep.q_elscales[i + s * ctx.n_qelscales];
            ctx.q_elscales->cpu_data_p[i + s * ctx.n_qelscales] = {row_int(row, 0), row_int(row, 1), row_double(row, 2)};
        }
    }

    size_t q_excl_count = 0;
    for (size_t s = 0; s < fep.q_exclpairs.size(); s++) q_excl_count += fep.q_exclpairs[s].size();
    ctx.n_qexclpairs = ctx.n_lambdas > 0 ? (int) q_excl_count / ctx.n_lambdas : 0;
    ctx.q_exclpairs.resize(ctx.n_qexclpairs * ctx.n_lambdas);
    for (int s = 0; s < ctx.n_lambdas && (size_t) s < fep.q_exclpairs.size(); s++) {
        for (int i = 0; i < ctx.n_qexclpairs && (size_t) i < fep.q_exclpairs[s].size(); i++) {
            std::vector<std::string> row = split_ws(fep.q_exclpairs[s][i]);
            ctx.q_exclpairs[i + s * ctx.n_qexclpairs] = {row_int(row, 0), row_int(row, 1), row_int(row, 2)};
        }
    }
}

static double q_randm(int &seed) {
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
    double r = ((irand / 10) * 10) / (double) m;
    if (r <= 0.0 || r > 1.0) r = 0.0;
    seed = irand;
    return r;
}

static double q_gauss(int &seed, double sd) {
    double total = 0.0;
    for (int i = 0; i < 12; i++) total += q_randm(seed);
    return (total - 6.0) * sd;
}

static std::vector<std::vector<std::string> > generate_velocities(const TopData &top, const MDInput &input) {
    std::map<std::string, std::string> mdv = input.keyed.count("md") ? input.keyed.find("md")->second : std::map<std::string, std::string>();
    int seed = atoi(value_or(mdv, "random-seed", value_or(mdv, "random_seed", "1")).c_str());
    double temperature = atof(value_or(mdv, "initial-temperature", value_or(mdv, "initial_temperature", value_or(mdv, "temperature", "0"))).c_str());
    if (seed <= 0) fatal("Fresh-start Q input requires a positive random_seed.");

    std::vector<std::vector<std::string> > velocities;
    double kT = Boltz * temperature;
    for (size_t i = 0; i < top.atypes.size(); i++) {
        int type_code = top.atypes[i].second;
        std::map<int, std::vector<std::string> >::const_iterator catype = top.catypes.find(type_code);
        if (catype == top.catypes.end() || catype->second.empty()) fatal("Topology atom type references missing mass.");
        double mass = atof(catype->second[0].c_str());
        double sd = sqrt(kT / mass);
        std::vector<std::string> row;
        for (int axis = 0; axis < 3; axis++) {
            char buffer[80];
            snprintf(buffer, sizeof(buffer), "%.15g", q_gauss(seed, sd));
            row.push_back(buffer);
        }
        velocities.push_back(row);
    }
    return velocities;
}

static int int_setting(const std::map<std::string, std::string> &values, const std::string &key, const std::string &fallback) {
    return atoi(value_or(values, key, fallback).c_str());
}

static std::string join_space(const std::vector<std::string> &fields) {
    std::string out;
    for (size_t i = 0; i < fields.size(); i++) {
        if (i > 0) out += " ";
        out += fields[i];
    }
    return out;
}

static std::string trajectory_atoms_setting(const MDInput &input) {
    std::map<std::string, std::vector<std::vector<std::string> > >::const_iterator it = input.lists.find("trajectory-atoms");
    if (it == input.lists.end() || it->second.empty()) return "";
    if (it->second.size() > 1) fatal("Native QGPU output supports only one trajectory_atoms mask line.");
    return join_space(it->second[0]);
}

static bool read_fortran_record(std::ifstream &in, std::vector<char> &payload) {
    int32_t nbytes = 0;
    in.read((char*) &nbytes, sizeof(nbytes));
    if (!in) return false;
    if (nbytes <= 0) fatal("Invalid Fortran restart record marker.");
    payload.assign(nbytes, 0);
    in.read(payload.data(), nbytes);
    int32_t end_nbytes = 0;
    in.read((char*) &end_nbytes, sizeof(end_nbytes));
    if (!in || end_nbytes != nbytes) fatal("Invalid or truncated Fortran restart record.");
    return true;
}

static std::vector<std::vector<std::string> > unpack_restart_vector(const std::vector<char> &payload) {
    if (payload.size() < sizeof(int32_t)) fatal("Invalid restart vector record.");
    int32_t nat3 = 0;
    memcpy(&nat3, payload.data(), sizeof(int32_t));
    if (nat3 <= 0 || nat3 % 3 != 0) fatal("Invalid restart atom count.");
    size_t expected = sizeof(int32_t) + (size_t) nat3 * sizeof(double);
    if (payload.size() != expected) fatal("Invalid restart vector record length.");
    const double *values = (const double*) (payload.data() + sizeof(int32_t));
    std::vector<std::vector<std::string> > rows;
    for (int i = 0; i < nat3; i += 3) {
        std::vector<std::string> row;
        for (int j = 0; j < 3; j++) {
            char buffer[80];
            snprintf(buffer, sizeof(buffer), "%.15g", values[i + j]);
            row.push_back(buffer);
        }
        rows.push_back(row);
    }
    return rows;
}

static void read_restart_csv_inputs(const std::string &path, std::vector<std::vector<std::string> > &coords, std::vector<std::vector<std::string> > &velocities) {
    std::ifstream in(path.c_str(), std::ios::binary);
    if (!in) fatal("Could not open restart file " + path);
    std::vector<char> payload;
    if (!read_fortran_record(in, payload)) fatal("Could not read restart coordinates.");
    coords = unpack_restart_vector(payload);
    if (!read_fortran_record(in, payload)) fatal("Could not read restart velocities.");
    velocities = unpack_restart_vector(payload);
    if (coords.size() != velocities.size()) fatal("Restart coordinate and velocity atom counts differ.");
}

void q_prepare_input(const char *input_file, const char *workdir) {
    std::string input_path(input_file);
    (void) workdir;
    std::string input_dir = dirname_of(input_path);

    MDInput input = parse_md_input(input_path);
    if (!input.keyed.count("files")) fatal("Q input requires a [files] section.");
    std::map<std::string, std::string> files = input.keyed.find("files")->second;
    std::string topology_file = resolve_existing_or_cwd(value_or(files, "topology", ""), input_dir);
    std::string fep_file = resolve_existing_or_cwd(value_or(files, "fep", ""), input_dir);
    std::string restart_file = resolve_existing_or_cwd(value_or(files, "restart", ""), input_dir);
    std::string final_file = resolve_output_path(value_or(files, "final", ""), input_dir);
    std::string trajectory_file = resolve_output_path(value_or(files, "trajectory", ""), input_dir);
    std::string energy_file = resolve_output_path(value_or(files, "energy", ""), input_dir);
    std::string trajectory_atoms = trajectory_atoms_setting(input);
    std::map<std::string, std::string> mdv = input.keyed.count("md") ? input.keyed.find("md")->second : std::map<std::string, std::string>();
    std::map<std::string, std::string> intervals = input.keyed.count("intervals") ? input.keyed.find("intervals")->second : std::map<std::string, std::string>();
    int steps = int_setting(mdv, "steps", "0");
    int trajectory_interval = int_setting(intervals, "trajectory", "0");
    int energy_interval = int_setting(intervals, "energy", "0");

    if (topology_file.empty()) fatal("Q input [files] section must specify topology.");
    q_output_configure(final_file.c_str(),
                       trajectory_file.c_str(),
                       energy_file.c_str(),
                       trajectory_atoms.c_str(),
                       topology_file.c_str(),
                       steps,
                       trajectory_interval,
                       energy_interval);
    TopData top = parse_topology(topology_file);
    FepData fep = parse_fep(fep_file);
    complete_fep_defaults(fep, top);

    std::vector<std::vector<std::string> > input_coords;
    std::vector<std::vector<std::string> > input_velocities;
    q_native_fresh_start = restart_file.empty();
    if (!restart_file.empty()) {
        read_restart_csv_inputs(restart_file, input_coords, input_velocities);
    }
    else {
        input_coords = top.coords;
        input_velocities = generate_velocities(top, input);
    }
    if (input_coords.size() != top.coords.size() || input_velocities.size() != top.coords.size()) {
        fatal("Input coordinate/velocity atom count does not match topology.");
    }

    set_md_context(input);
    set_topology_context(top, input_coords);
    set_fep_context(fep);
    set_velocities_context(input_velocities);
}
