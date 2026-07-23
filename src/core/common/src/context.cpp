#include "context.h"

#include <algorithm>
#include <cstdio>
#include <cstdlib>

#include "constants.h"
#include "geometry.h"
#include "helpers.h"
#include "inp_parser.h"
#include "parse.h"
#include "shake.h"

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

template <typename T>
bool same_size(
    const std::vector<T>& a,
    const std::vector<T>& b) {
    return a.size() == b.size();
}

bool same_coords(
    const std::vector<coord_t>& a,
    const std::vector<coord_t>& b) {
    if (a.size() != b.size()) {
        return false;
    }

    for (size_t i = 0; i < a.size(); i++) {
        if (a[i].x != b[i].x ||
            a[i].y != b[i].y ||
            a[i].z != b[i].z) {
            return false;
        }
    }

    return true;
}

bool same_topology_scalars(
    const topo_t& a,
    const topo_t& b) {
    return a.solvent_type == b.solvent_type &&
           a.exclusion_radius == b.exclusion_radius &&
           a.solvent_radius == b.solvent_radius &&

           a.solute_center.x == b.solute_center.x &&
           a.solute_center.y == b.solute_center.y &&
           a.solute_center.z == b.solute_center.z &&

           a.solvent_center.x == b.solvent_center.x &&
           a.solvent_center.y == b.solvent_center.y &&
           a.solvent_center.z == b.solvent_center.z &&

           a.el14_scale == b.el14_scale &&
           a.coulomb_constant == b.coulomb_constant &&
           a.vdw_rule == b.vdw_rule;
}

bool same_md_batch_config(
    const md_t& a,
    const md_t& b) {
    return
        // Integration and thermostat
        a.steps == b.steps &&
        a.stepsize == b.stepsize &&
        a.temperature == b.temperature &&
        a.thermostat == b.thermostat &&
        a.bath_coupling == b.bath_coupling &&
        a.initial_temperature == b.initial_temperature &&

        // SHAKE
        a.shake_solvent == b.shake_solvent &&
        a.shake_solute == b.shake_solute &&
        a.shake_hydrogens == b.shake_hydrogens &&

        // Nonbonded configuration
        a.lrf == b.lrf &&
        a.separate_scaling == b.separate_scaling &&
        a.charge_groups == b.charge_groups &&
        a.solute_solute == b.solute_solute &&
        a.solvent_solvent == b.solvent_solvent &&
        a.solute_solvent == b.solute_solvent &&
        a.q_atom == b.q_atom &&

        // Spherical boundary
        a.shell_radius == b.shell_radius &&
        a.shell_force == b.shell_force &&
        a.radial_force == b.radial_force &&
        a.polarisation == b.polarisation &&
        a.charge_correction == b.charge_correction &&
        a.polarisation_force == b.polarisation_force &&

        // Runtime/output intervals
        a.non_bond == b.non_bond &&
        a.output == b.output &&
        a.energy == b.energy &&
        a.trajectory == b.trajectory;
}

void validate_batch_compatible(
    const ParseResult& reference,
    const ParseResult& candidate,
    int replica) {
    const bool same_run_configuration =
        reference.fresh_start == candidate.fresh_start &&
        same_md_batch_config(reference.md, candidate.md) &&
        same_topology_scalars(reference.topo, candidate.topo);

    const bool same_system_shape =
        reference.n_atoms_solute ==
            candidate.n_atoms_solute &&
        reference.n_bonds_solute ==
            candidate.n_bonds_solute &&
        reference.n_angles_solute ==
            candidate.n_angles_solute &&
        reference.n_torsions_solute ==
            candidate.n_torsions_solute &&
        reference.n_impropers_solute ==
            candidate.n_impropers_solute &&
        reference.lambdas == candidate.lambdas &&
        same_coords(
            reference.coords_init,
            candidate.coords_init);

    const bool same_classical_topology =
        same_size(reference.bonds, candidate.bonds) &&
        same_size(reference.cbonds, candidate.cbonds) &&
        same_size(reference.angles, candidate.angles) &&
        same_size(reference.cangles, candidate.cangles) &&
        same_size(reference.torsions, candidate.torsions) &&
        same_size(reference.ctorsions, candidate.ctorsions) &&
        same_size(reference.impropers, candidate.impropers) &&
        same_size(
            reference.cimpropers,
            candidate.cimpropers) &&
        same_size(reference.charges, candidate.charges) &&
        same_size(reference.ccharges, candidate.ccharges) &&
        same_size(reference.atypes, candidate.atypes) &&
        same_size(reference.catypes, candidate.catypes) &&
        reference.excluded == candidate.excluded &&
        reference.molecules == candidate.molecules &&
        reference.ngbrs14 == candidate.ngbrs14 &&
        reference.ngbrs14_long == candidate.ngbrs14_long &&
        reference.ngbrs23 == candidate.ngbrs23 &&
        reference.ngbrs23_long == candidate.ngbrs23_long;

    const bool same_restraints =
        same_size(
            reference.restrspos,
            candidate.restrspos) &&
        same_size(
            reference.restrangs,
            candidate.restrangs) &&
        same_size(
            reference.restrdists,
            candidate.restrdists) &&
        same_size(
            reference.restrseqs,
            candidate.restrseqs) &&
        same_size(
            reference.restrwalls,
            candidate.restrwalls);

    const bool same_fep_topology =
        reference.q_atoms == candidate.q_atoms &&
        same_size(
            reference.q_angcouples,
            candidate.q_angcouples) &&
        same_size(reference.q_atypes, candidate.q_atypes) &&
        same_size(reference.q_catypes, candidate.q_catypes) &&
        same_size(reference.q_charges, candidate.q_charges) &&
        same_size(reference.q_angles, candidate.q_angles) &&
        same_size(reference.q_cangles, candidate.q_cangles) &&
        same_size(reference.q_bonds, candidate.q_bonds) &&
        same_size(reference.q_cbonds, candidate.q_cbonds) &&
        same_size(
            reference.q_impropers,
            candidate.q_impropers) &&
        same_size(
            reference.q_cimpropers,
            candidate.q_cimpropers) &&
        same_size(
            reference.q_torsions,
            candidate.q_torsions) &&
        same_size(
            reference.q_ctorsions,
            candidate.q_ctorsions) &&
        same_size(
            reference.q_offdiags,
            candidate.q_offdiags) &&
        same_size(
            reference.q_imprcouples,
            candidate.q_imprcouples) &&
        same_size(
            reference.q_softpairs,
            candidate.q_softpairs) &&
        same_size(
            reference.q_torcouples,
            candidate.q_torcouples) &&
        same_size(
            reference.q_elscales,
            candidate.q_elscales) &&
        same_size(
            reference.q_exclpairs,
            candidate.q_exclpairs) &&
        same_size(
            reference.q_shakes,
            candidate.q_shakes) &&
        same_size(
            reference.q_softcores,
            candidate.q_softcores) &&
        reference.softcore_use_max_potential ==
            candidate.softcore_use_max_potential;

    // Coordinate and velocity values may differ, but every replica
    // must contain the same number of atoms.
    const bool same_dynamic_state_shape =
        reference.coords.size() ==
            candidate.coords.size() &&
        reference.velocities.size() ==
            candidate.velocities.size();

    const bool compatible =
        same_run_configuration &&
        same_system_shape &&
        same_classical_topology &&
        same_restraints &&
        same_fep_topology &&
        same_dynamic_state_shape;

    if (!compatible) {
        fatal(
            "Replica " + std::to_string(replica) +
            " is not topology/FEP compatible with replica 0.");
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

std::vector<ParseResult> Context::get_parse_results() {
    const auto& input_files = command_info.input_files;
    if (input_files.empty()) {
        fatal("At least one replica input file is required.");
    }

    std::vector<ParseResult> results;
    results.reserve(input_files.size());

    for (const std::string& input_file : command_info.input_files) {
        InpParser replica_parser(input_file);
        ParseResult candidate = replica_parser.parse();
        validate_parse_result(candidate);
        const int replica = results.size();
        if (!results.empty()) {
            validate_batch_compatible(results.front(), candidate, replica);
        }
        results.push_back(std::move(candidate));
    }
    return results;
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

void Context::init_native_outputs(const std::vector<ParseResult>& replica_results) {
    native_outputs.clear();
    native_outputs.reserve(replica_results.size());

    for (const ParseResult& result : replica_results) {
        native_outputs.push_back(result.native_output);
    }
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

void Context::init_coords(
    const std::vector<ParseResult>& replica_results) {
    std::vector<coord_t> batch_coords;

    const size_t parsed_atoms = static_cast<size_t>(n_atoms) * replica_results.size();

    batch_coords.reserve(parsed_atoms);

    for (const ParseResult& result : replica_results) {
        batch_coords.insert(batch_coords.end(), result.coords.begin(), result.coords.end());
    }

    coords = HostDeviceBuffer<coord_t>::from_vector(batch_coords, command_info.requested_gpu);
}

void Context::init_velocities(const std::vector<ParseResult>& replica_results) {
    std::vector<vel_t> batch_velocities;

    batch_velocities.reserve(static_cast<size_t>(n_atoms) * replica_results.size());

    for (const ParseResult& result : replica_results) {
        batch_velocities.insert(batch_velocities.end(), result.velocities.begin(), result.velocities.end());
    }

    velocities = HostDeviceBuffer<vel_t>::from_vector(batch_velocities, command_info.requested_gpu);
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
    dvelocities = std::make_unique<HostDeviceBuffer<dvel_t>>(n_total_atoms(), true, command_info.requested_gpu);
}

void Context::init_energy() {
    energy.init(n_lambdas(), n_replicates());
}

void Context::init(const std::vector<ParseResult>& replica_results) {
    if (replica_results.empty()) {
        fatal("At least one parsed replica is required.");
    }
    const ParseResult& parsed = replica_results.front();

    // 1. Basic configuration
    init_fresh_start(parsed);
    init_md(parsed);
    init_topo(parsed);
    init_native_outputs(replica_results);
    init_const(parsed);

    // 2. Basic atom data
    init_coords_init(parsed);
    init_coords(replica_results);
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
    init_velocities(replica_results);

    // 6. Runtime buffers
    init_dvelocities();
    init_energy();
}