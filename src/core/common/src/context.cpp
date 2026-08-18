#include "context.h"

#include <algorithm>
#include <cstdio>
#include <cstdlib>

#include "constants.h"
#include "geometry.h"
#include "helpers.h"
#include "inp_parser.h"
#include "parse.h"
#include "q_math.h"

namespace {

int per_lambda_count(size_t size, int n_lambdas) {
    return n_lambdas > 0 ? static_cast<int>(size) / n_lambdas : 0;
}

void validate_parse_result(const ParseResult& parsed) {
    const size_t n_atoms = parsed.coords_init.size();
    if (parsed.coords.size() != n_atoms || parsed.velocities.size() != n_atoms) {
        fatal("Parsed coordinate, initial coordinate, and velocity counts must match.");
    }
    if (parsed.excluded.size() != n_atoms) {
        fatal("Parsed excluded flags must match atom count.");
    }
    if (parsed.topo.vdw_rule != VDW_GEOMETRIC && parsed.topo.vdw_rule != VDW_ARITHMETIC) {
        fatal("Invalid vdw_rule in parsed topology.");
    }
}

}  // namespace

void Context::set_lj_pairs(std::vector<int>& matrix, const std::vector<std::pair<int, int>>& pairs, int value) {
    for (const auto& pair : pairs) {
        int ai = pair.first;
        int aj = pair.second;
        if (ai >= 0 && aj >= 0 && ai < n_atoms_solute && aj < n_atoms_solute) {
            matrix[ai * n_atoms_solute + aj] = value;
            matrix[aj * n_atoms_solute + ai] = value;
        }
    }
}

void Context::init_lj_matrix(const ParseResult& parsed) {
    std::vector<int> matrix(n_atoms_solute * n_atoms_solute, 0);
    set_lj_pairs(matrix, parsed.ngbrs14, 1);
    set_lj_pairs(matrix, parsed.ngbrs14_long, 1);
    set_lj_pairs(matrix, parsed.ngbrs23, 3);
    set_lj_pairs(matrix, parsed.ngbrs23_long, 3);

    LJ_matrix = HostDeviceBuffer<int>::from_vector(matrix, command_info.requested_gpu);
}

ParseResult Context::get_parse_result() {
    std::unique_ptr<BaseParser> parser = std::make_unique<InpParser>(command_info.input_file);
    ParseResult parser_result = parser->parse();
    return parser_result;
}

void Context::init_fresh_start(const ParseResult& parsed) {
    fresh_start = parsed.fresh_start;
}

void Context::init_md(const ParseResult& parsed) {
    md = parsed.md;
}

void Context::init_topo(const ParseResult& parsed) {
    topo = parsed.topo;
}

void Context::init_native_output(const ParseResult& parsed) {
    native_output = parsed.native_output;
}

void Context::init_const(const ParseResult& parsed) {
    // Need data from md
    n_atoms = static_cast<int>(parsed.coords_init.size());
    n_atoms_solute = parsed.n_atoms_solute;

    dt = time_unit * md.stepsize;
}

void Context::init_coords_init(const ParseResult& parsed) {
    coords_init = HostDeviceBuffer<coord_t>::from_vector(parsed.coords_init, command_info.requested_gpu);
}

void Context::init_coords(const ParseResult& parsed) {
    coords = HostDeviceBuffer<coord_t>::from_vector(parsed.coords, command_info.requested_gpu);
}

void Context::init_velocities(const ParseResult& parsed) {
    if (!fresh_start) {
        velocities = HostDeviceBuffer<vel_t>::from_vector(parsed.velocities, command_info.requested_gpu);
    } else {
        srand(md.random_seed);
        const auto* atom_types = atypes->cpu_data_p;
        const auto* type_parameters = catypes->cpu_data_p;

        std::vector<vel_t> velocity_data(n_atoms);

        const double kT = Boltz * md.initial_temperature;
        for (int atom = 0; atom < n_atoms; atom++) {
            double mass = type_parameters[atom_types[atom].code - 1].m;
            const double sd = sqrt(kT / mass);
            velocity_data[atom].x = gauss(0, sd);
            velocity_data[atom].y = gauss(0, sd);
            velocity_data[atom].z = gauss(0, sd);
        }
        velocities = HostDeviceBuffer<vel_t>::from_vector(velocity_data, command_info.requested_gpu);
    }
}

void Context::init_lambdas(const ParseResult& parsed) {
    lambdas = HostDeviceBuffer<double>::from_vector(parsed.lambdas, command_info.requested_gpu);
}

void Context::init_charges(const ParseResult& parsed) {
    charges = HostDeviceBuffer<charge_t>::from_vector(parsed.charges, command_info.requested_gpu);
}

void Context::init_ccharges(const ParseResult& parsed) {
    ccharges = HostDeviceBuffer<ccharge_t>::from_vector(parsed.ccharges, command_info.requested_gpu);
}

void Context::init_atypes(const ParseResult& parsed) {
    atypes = HostDeviceBuffer<atype_t>::from_vector(parsed.atypes, command_info.requested_gpu);
}

void Context::init_catypes(const ParseResult& parsed) {
    if (topo.vdw_rule == VDW_GEOMETRIC) {
        catypes = HostDeviceBuffer<catype_t>::from_vector(parsed.catypes, command_info.requested_gpu);
    } else {
        auto catypes_ = parsed.catypes;
        int n_catypes = catypes_.size();
        for (int i = 0; i < n_catypes; i++) {
            catype_t& catype = catypes_[i];
            catype.bii_normal = sqrt(abs(catype.bii_normal));
            catype.bii_1_4 = sqrt(abs(catype.bii_1_4));
        }
        catypes = HostDeviceBuffer<catype_t>::from_vector(catypes_, command_info.requested_gpu);
    }
}

void Context::init_excluded(const ParseResult& parsed) {
    excluded = HostDeviceBuffer<bool>::from_vector(parsed.excluded, command_info.requested_gpu);

    n_excluded = 0;
    for (const bool is_excluded : parsed.excluded) {
        if (is_excluded) {
            ++n_excluded;
        }
    }
}

void Context::init_molecules(const ParseResult& parsed) {
    molecules = parsed.molecules;
}

void Context::init_charge_group_config(const ParseResult& parsed) {
    charge_group_config = parsed.charge_groups;
}

void Context::init_q_atoms(const ParseResult& parsed) {
    q_atoms = parsed.q_atoms;
}

void Context::init_q_catypes(const ParseResult& parsed) {
    q_catypes = parsed.q_catypes;
    if (topo.vdw_rule == VDW_GEOMETRIC) {
    } else {
        for (catype_t& catype : q_catypes) {
            catype.bii_normal = sqrt(abs(catype.bii_normal));
            catype.bii_1_4 = sqrt(abs(catype.bii_1_4));
        }
    }
}

void Context::init_q_atypes(const ParseResult& parsed) {
    q_atypes = parsed.q_atypes;
}

void Context::init_q_charges(const ParseResult& parsed) {
    q_charges = parsed.q_charges;
}

void Context::init_inv_mass() {
    const auto* atom_types = atypes->cpu_data_p;
    const auto* type_parameters = catypes->cpu_data_p;
    std::vector<double> inv_w(n_atoms);
    for (int atom = 0; atom < n_atoms; atom++) {
        inv_w[atom] = 1.0 / type_parameters[atom_types[atom].code - 1].m;
    }
    winv = HostDeviceBuffer<double>::from_vector(inv_w, command_info.requested_gpu);
}

void Context::init_heavy() {
    const auto& host_catypes = catypes->cpu_data_p;
    const auto& host_atypes = atypes->cpu_data_p;
    std::vector<bool> heavy_(n_atoms, false);
    for (int i = 0; i < n_atoms; i++) {
        double mass = host_catypes[host_atypes[i].code - 1].m;
        heavy_[i] = (mass >= 4.0);
    }
    heavy = HostDeviceBuffer<bool>::from_vector(heavy_, command_info.requested_gpu);
}

void Context::init_patoms() {
    int n_patoms = n_atoms_solute - n_qatoms();
    p_atoms.resize(n_patoms);

    int pi = 0, qi = 0;
    for (int i = 0; i < n_atoms_solute; i++) {
        if (qi < n_qatoms() && i == q_atoms[qi]) {
            qi++;
        } else {
            p_atoms[pi] = i;
            pi++;
        }
    }
}

void Context::init_shell() {
    const auto* host_heavy = heavy->cpu_data_p;
    const auto* host_excluded = excluded->cpu_data_p;
    const auto* host_coords_init = coords_init->cpu_data_p;

    std::vector<bool> shell_(n_atoms, false);

    if (md.charge_groups) {
        double rin2 = md.shell_radius * md.shell_radius;
        auto& host_charge_groups = charge_group_config;
        const bool use_switch_atom = host_charge_groups.iuse_switch_atom == 1;

        for (int grp = 0; grp < host_charge_groups.n_cgrps_solute; grp++) {
            const auto& charge_group = host_charge_groups.charge_groups[grp];
            int i = charge_group.iswitch - 1;
            if (i >= 0 && i < n_atoms_solute && host_heavy[i] && !host_excluded[i]) {
                coord_t c = host_coords_init[i];
                if (!use_switch_atom) {
                    c = {0, 0, 0};
                    for (int atom : charge_group.atoms) {
                        int ai = atom - 1;
                        c = c + host_coords_init[ai];
                    }
                    double inv_atoms = 1.0 / static_cast<double>(charge_group.atoms.size());
                    c = c * inv_atoms;
                }

                double r2 = norm2(c - topo.solute_center);
                bool in_shell = r2 > rin2;
                for (int atom : charge_group.atoms) {
                    shell_[atom - 1] = in_shell;
                }
            }
        }
    } else {
        double rin2 = pow(shell_default * topo.exclusion_radius, 2);

        for (int i = 0; i < n_atoms; i++) {
            if (i >= 0 && i < n_atoms_solute && host_heavy[i] && !host_excluded[i]) {
                double dist2 = norm2(host_coords_init[i] - topo.solute_center);
                if (dist2 > rin2) {
                    shell_[i] = true;
                }
            }
        }
    }
    shell = HostDeviceBuffer<bool>::from_vector(shell_, command_info.requested_gpu);
}

void Context::init_dvelocities() {
    dvelocities = std::make_unique<HostDeviceBuffer<dvel_t>>(n_atoms, true, command_info.requested_gpu);
}

void Context::init_energy() {
    energy.init(n_lambdas());
}

void Context::init(const ParseResult& parsed) {
    validate_parse_result(parsed);

    // 1. Basic configuration
    init_fresh_start(parsed);
    init_md(parsed);
    init_topo(parsed);
    init_native_output(parsed);
    init_const(parsed);

    // 2. Basic atom data
    init_coords_init(parsed);
    init_coords(parsed);
    init_charges(parsed);
    init_ccharges(parsed);
    init_atypes(parsed);
    init_catypes(parsed);
    init_excluded(parsed);
    init_molecules(parsed);
    init_charge_group_config(parsed);

    // 3. FEP data
    init_lambdas(parsed);
    init_q_atoms(parsed);
    init_q_atypes(parsed);
    init_q_catypes(parsed);
    init_q_charges(parsed);

    // 4. Derived context data
    init_inv_mass();
    init_heavy();
    init_shell();
    init_patoms();
    init_lj_matrix(parsed);

    // 5. Dynamic state
    init_velocities(parsed);

    // 6. Runtime buffers
    init_dvelocities();
    init_energy();
}