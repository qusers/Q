#include "context.h"

#include <algorithm>
#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <string>

#include "constants.h"
#include "geometry.h"
#include "helpers.h"
#include "inp_parser.h"
#include "parse.h"
#include "q_math.h"
#include "qgpu_error.h"
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
        velocities = std::make_unique<HostDeviceBuffer<vel_t>>(n_atoms, true, command_info.requested_gpu);
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

void Context::init_offdiags(const ParseResult& parsed) {
    const bool run_gpu = command_info.requested_gpu;

    std::vector<offdiag_term_t> terms;
    std::vector<offdiag_record_t> records;

    terms.reserve(parsed.q_offdiags.size());
    records.reserve(parsed.q_offdiags.size());

    for (const auto& param : parsed.q_offdiags) {
        // q_atoms already contains 0-based physical atom indices.
        const int ak = q_atoms.at(param.qk - 1);
        const int al = q_atoms.at(param.ql - 1);

        offdiag_term_t term{};
        term.state_i = param.i;
        term.state_j = param.j;
        term.ak = ak;
        term.al = al;
        term.A = param.Aij;
        term.mu = param.muij;
        terms.push_back(term);

        offdiag_record_t record{};
        record.state_i = param.i;
        record.state_j = param.j;
        records.push_back(record);
    }

    offdiag_terms = HostDeviceBuffer<offdiag_term_t>::from_vector(terms, run_gpu);

    offdiag_records = HostDeviceBuffer<offdiag_record_t>::from_vector(records, run_gpu);
}

void Context::init_masses(const ParseResult& parsed) {
    const auto* atom_types = atypes->cpu_data_p;
    const auto* type_parameters = catypes->cpu_data_p;
    const auto* atom_heavy = heavy->cpu_data_p;

    std::vector<double> atom_masses(n_atoms);
    for (int i = 0; i < n_atoms; i++) {
        const int type = atom_types[i].code - 1;
        const double mass = type_parameters[type].m;

        atom_masses[i] = mass;
    }

    if (md.hmr) {
        std::vector<int> hydrogen_parent(n_atoms, -1);
        const double target = md.hmr_target_mass;
        int repartitioned_hydrogens = 0;

        for (const bond_t& bond : parsed.bonds) {
            const int i = bond.ai - 1;
            const int j = bond.aj - 1;

            // Solute only HMR: ignore solvent bond
            if (i >= n_atoms_solute || j >= n_atoms_solute) {
                continue;
            }
            const bool ih = atom_heavy[i];
            const bool jh = atom_heavy[j];
            if (ih == jh) {
                continue;
            }

            const int hydrogen = ih ? j : i;
            const int parent = ih ? i : j;

            if (hydrogen_parent[hydrogen] != -1 && hydrogen_parent[hydrogen] != parent) {
                fatal("HMR hydrogen " + std::to_string(hydrogen + 1) + " is bonded to multiple heavy atoms");
            }
            hydrogen_parent[hydrogen] = parent;
            const double old_mass = atom_masses[hydrogen];
            const double delta = old_mass - target;

            atom_masses[hydrogen] = target;
            atom_masses[parent] += delta;
            repartitioned_hydrogens++;
        }
        printf("HMR enabled! Target hydrogen mass: 3.024 amu. Repartitioned solute hydrogens: %d\n", repartitioned_hydrogens);
    }
    masses = HostDeviceBuffer<double>::from_vector(atom_masses, command_info.requested_gpu);
}

void Context::init_inv_mass() {
    const double* atom_masses = masses->cpu_data_p;
    std::vector<double> inv_w(n_atoms);
    for (int atom = 0; atom < n_atoms; atom++) {
        inv_w[atom] = 1.0 / atom_masses[atom];
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
    // q_atoms stores 0-based physical atom indices.
    std::vector<bool> is_q(n_atoms_solute, false);
    for (const int atom : q_atoms) {
        is_q.at(atom) = true;
    }

    p_atoms.clear();
    p_atoms.reserve(n_atoms_solute - n_qatoms());

    for (int atom = 0; atom < n_atoms_solute; ++atom) {
        if (!is_q[atom]) {
            p_atoms.push_back(atom);
        }
    }
}

void Context::init_shell() {
    const auto* host_heavy = heavy->cpu_data_p;
    const auto* host_excluded = excluded->cpu_data_p;
    const auto* host_coords_init = coords_init->cpu_data_p;

    std::vector<bool> shell_(n_atoms, false);

    if (md.charge_groups) {
        // As in Fortran Q: shell_radius <= 1 is a fraction of the exclusion radius, else Å.
        if (md.shell_radius > topo.exclusion_radius) {
            throw InputError("shell_radius must be < the exclusion radius (" + std::to_string(topo.exclusion_radius) + ").");
        }
        const double rin = md.shell_radius <= 1.0 ? md.shell_radius * topo.exclusion_radius : md.shell_radius;
        double rin2 = rin * rin;
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
    init_offdiags(parsed);

    // 4. Derived context data
    init_heavy();
    init_masses(parsed);
    init_inv_mass();
    init_shell();
    init_patoms();
    init_lj_matrix(parsed);

    // 5. Dynamic state
    init_velocities(parsed);

    // 6. Runtime buffers
    init_dvelocities();
    init_energy();
}

void Context::generate_velocities(int seed, double initial_temperature) {
    md.random_seed = seed;
    md.initial_temperature = initial_temperature;
    srand(seed);
    const auto* atom_mass = masses->cpu_data_p;

    std::vector<vel_t> velocity_data(n_atoms);
    const double kT = Boltz * initial_temperature;
    for (int atom = 0; atom < n_atoms; atom++) {
        const double sd = sqrt(kT / atom_mass[atom]);
        velocity_data[atom].x = gauss(0, sd);
        velocity_data[atom].y = gauss(0, sd);
        velocity_data[atom].z = gauss(0, sd);
    }

    if (!velocities || velocities->length != static_cast<size_t>(n_atoms)) {
        throw InputError("generate_velocities: velocity buffer not allocated for " + std::to_string(n_atoms) + " atoms.");
    }
    std::copy(velocity_data.begin(), velocity_data.end(), velocities->cpu_data_p);
    if (command_info.requested_gpu) velocities->upload();
    fresh_start = true;
}

void Context::validate_system_fields(const md_t& s) const {
    auto check = [](bool same, const char* name) {
        if (!same) {
            throw InputError(std::string("Stage changes system-level setting '") + name +
                             "'; it must stay the same for the whole engine.");
        }
    };
    check(s.hmr == md.hmr, "hmr");
    check(s.hmr_target_mass == md.hmr_target_mass, "hmr_target_mass");
    check(s.shell_radius == md.shell_radius, "shell_radius");
    check(s.charge_groups == md.charge_groups, "charge_groups");
    check(s.lrf == md.lrf, "lrf");
    check(s.lrf_cutoff == md.lrf_cutoff, "lrf_cutoff");
    check(s.solute_solute == md.solute_solute, "solute_solute");
    check(s.solvent_solvent == md.solvent_solvent, "solvent_solvent");
    check(s.solute_solvent == md.solute_solvent, "solute_solvent");
    check(s.q_atom == md.q_atom, "q_atom");
    check(s.polarisation == md.polarisation, "polarisation");
    check(s.charge_correction == md.charge_correction, "charge_correction");
    check(s.radial_force == md.radial_force, "radial_force");
    check(s.polarisation_force == md.polarisation_force, "polarisation_force");
    check(s.exclude_bonded == md.exclude_bonded, "exclude_bonded");
    check(s.non_bond == md.non_bond, "non_bond");
    check(s.solute_constraint_algorithm == md.solute_constraint_algorithm, "solute_constraint_algorithm");
    check(s.solvent_constraint_algorithm == md.solvent_constraint_algorithm, "solvent_constraint_algorithm");
}

void Context::validate_stage(const StageSettings& stage) const {
    const md_t& s = stage.md;
    validate_system_fields(s);
    if (stage.lambdas.size() != static_cast<size_t>(n_lambdas())) {
        throw InputError("Stage has " + std::to_string(stage.lambdas.size()) +
                         " lambdas, system has " + std::to_string(n_lambdas()) + " FEP states.");
    }
    if (s.steps < 0) throw InputError("Stage steps must be >= 0.");
    if (s.stepsize <= 0) throw InputError("Stage stepsize must be > 0.");
    if (s.temperature <= 0) throw InputError("Stage temperature must be > 0 (no dynamics at zero temperature).");
    if (s.bath_coupling < s.stepsize) throw InputError("Stage bath_coupling must be >= stepsize.");
}

void Context::apply_stage(const StageSettings& stage) {
    validate_stage(stage);  // before anything changes, so a rejected stage leaves the context as it was
    const md_t& s = stage.md;

    // Stage fields only; System fields in md stay as loaded.
    md.steps = s.steps;
    md.stepsize = s.stepsize;
    md.temperature = s.temperature;
    md.thermostat = s.thermostat;
    md.bath_coupling = s.bath_coupling;
    md.constraint_solvent = s.constraint_solvent;
    md.constraint_solute = s.constraint_solute;
    md.constraint_hydrogens = s.constraint_hydrogens;
    md.shake_solvent = s.constraint_solvent;  // legacy mirrors
    md.shake_solute = s.constraint_solute;
    md.shake_hydrogens = s.constraint_hydrogens;
    md.separate_scaling = s.separate_scaling;
    md.shell_force = s.shell_force;
    md.output = s.output;
    md.energy = s.energy;
    md.trajectory = s.trajectory;

    dt = time_unit * md.stepsize;

    // Lambdas: overwrite values in place; length is fixed by the FEP states.
    std::copy(stage.lambdas.begin(), stage.lambdas.end(), lambdas->cpu_data_p);
    if (command_info.requested_gpu) lambdas->upload();

    // Output paths for this stage; keep topology_file from the system.
    const std::string topology_file = native_output.topology_file;
    native_output = stage.outputs;
    native_output.topology_file = topology_file;
}

void Context::set_state(const SystemState& state) {
    if (state.coords.size() != static_cast<size_t>(n_atoms) ||
        state.velocities.size() != static_cast<size_t>(n_atoms)) {
        throw InputError("State atom count does not match the system (" + std::to_string(n_atoms) + ").");
    }
    std::copy(state.coords.begin(), state.coords.end(), coords->cpu_data_p);
    std::copy(state.velocities.begin(), state.velocities.end(), velocities->cpu_data_p);
    if (command_info.requested_gpu) {
        coords->upload();
        velocities->upload();
    }
    fresh_start = false;
}