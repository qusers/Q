#include "context.h"

#include <algorithm>
#include <cstdio>
#include <cstdlib>

#include "constants.h"
#include "csv_parser.h"
#include "helpers.h"
#include "init.h"
#include "inp_parser.h"
#include "parse.h"

Context* Context::current_ = nullptr;

namespace {
template <typename T>
std::unique_ptr<HostDeviceBuffer<T>> buffer_from_vector(const std::vector<T>& values, bool run_gpu) {
    auto buffer = std::make_unique<HostDeviceBuffer<T>>(values.size(), true, run_gpu);
    if (!values.empty()) {
        std::copy(values.begin(), values.end(), buffer->cpu_data_p);
    }
    if (run_gpu) {
        buffer->upload();
    }
    return buffer;
}

std::unique_ptr<HostDeviceBuffer<bool>> bool_buffer_from_vector(const std::vector<bool>& values, bool run_gpu) {
    auto buffer = std::make_unique<HostDeviceBuffer<bool>>(values.size(), true, run_gpu);
    for (size_t i = 0; i < values.size(); i++) {
        buffer->cpu_data_p[i] = values[i];
    }
    if (run_gpu) {
        buffer->upload();
    }
    return buffer;
}

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

void set_lj_pairs(Context& ctx, const std::vector<std::pair<int, int>>& pairs, int value) {
    int* matrix = ctx.LJ_matrix->cpu_data_p;
    for (const auto& pair : pairs) {
        int ai = pair.first;
        int aj = pair.second;
        if (ai >= 0 && aj >= 0 && ai < ctx.n_atoms_solute && aj < ctx.n_atoms_solute) {
            matrix[ai * ctx.n_atoms_solute + aj] = value;
            matrix[aj * ctx.n_atoms_solute + ai] = value;
        }
    }
}

void set_lj_matrix(Context& ctx, const ParseResult& parsed) {
    ctx.LJ_matrix = std::make_unique<HostDeviceBuffer<int>>(ctx.n_atoms_solute * ctx.n_atoms_solute, true, ctx.command_info.requested_gpu);
    set_lj_pairs(ctx, parsed.ngbrs14, 1);
    set_lj_pairs(ctx, parsed.ngbrs14_long, 1);
    set_lj_pairs(ctx, parsed.ngbrs23, 3);
    set_lj_pairs(ctx, parsed.ngbrs23_long, 3);
}

void preprocess_vdw_rule_parameters(Context& ctx) {
    if (ctx.topo.vdw_rule != VDW_ARITHMETIC) return;

    for (int i = 0; i < ctx.n_catypes; i++) {
        catype_t& catype = ctx.catypes->cpu_data_p[i];
        catype.bii_normal = sqrt(fabs(catype.bii_normal));
        catype.bii_1_4 = sqrt(fabs(catype.bii_1_4));
    }

    for (catype_t& catype : ctx.q_catypes) {
        catype.bii_normal = sqrt(fabs(catype.bii_normal));
        catype.bii_1_4 = sqrt(fabs(catype.bii_1_4));
    }

    if (ctx.command_info.requested_gpu) {
        ctx.catypes->upload();
    }
}

void upload_preprocessed_topology(Context& ctx) {
    if (!ctx.command_info.requested_gpu) return;

    ctx.bonds->upload();
    ctx.angles->upload();
    ctx.torsions->upload();
    ctx.impropers->upload();
}

void apply_parse_result(Context& ctx, const ParseResult& parsed) {
    validate_parse_result(parsed);

    ctx.fresh_start = parsed.fresh_start;
    bool run_gpu = ctx.command_info.requested_gpu;

    ctx.md = parsed.md;
    ctx.topo = parsed.topo;
    ctx.native_output = parsed.native_output;
    ctx.has_restart_wshell_theta_corr = parsed.has_restart_wshell_theta_corr;
    ctx.restart_wshell_theta_corr = parsed.restart_wshell_theta_corr;

    ctx.n_atoms = static_cast<int>(parsed.coords_init.size());
    ctx.n_atoms_solute = parsed.n_atoms_solute;
    ctx.n_bonds_solute = parsed.n_bonds_solute;
    ctx.n_angles_solute = parsed.n_angles_solute;
    ctx.n_torsions_solute = parsed.n_torsions_solute;
    ctx.n_impropers_solute = parsed.n_impropers_solute;

    ctx.coords_init = buffer_from_vector(parsed.coords_init, run_gpu);
    ctx.coords = buffer_from_vector(parsed.coords, run_gpu);
    ctx.velocities = buffer_from_vector(parsed.velocities, run_gpu);

    ctx.n_lambdas = static_cast<int>(parsed.lambdas.size());
    ctx.lambdas = buffer_from_vector(parsed.lambdas, run_gpu);

    ctx.n_bonds = static_cast<int>(parsed.bonds.size());
    ctx.bonds = buffer_from_vector(parsed.bonds, run_gpu);
    ctx.n_cbonds = static_cast<int>(parsed.cbonds.size());
    ctx.cbonds = buffer_from_vector(parsed.cbonds, run_gpu);

    ctx.n_angles = static_cast<int>(parsed.angles.size());
    ctx.angles = buffer_from_vector(parsed.angles, run_gpu);
    ctx.n_cangles = static_cast<int>(parsed.cangles.size());
    ctx.cangles = buffer_from_vector(parsed.cangles, run_gpu);

    ctx.n_torsions = static_cast<int>(parsed.torsions.size());
    ctx.torsions = buffer_from_vector(parsed.torsions, run_gpu);
    ctx.n_ctorsions = static_cast<int>(parsed.ctorsions.size());
    ctx.ctorsions = buffer_from_vector(parsed.ctorsions, run_gpu);

    ctx.n_impropers = static_cast<int>(parsed.impropers.size());
    ctx.impropers = buffer_from_vector(parsed.impropers, run_gpu);
    ctx.n_cimpropers = static_cast<int>(parsed.cimpropers.size());
    ctx.cimpropers = buffer_from_vector(parsed.cimpropers, run_gpu);

    ctx.n_restrspos = static_cast<int>(parsed.restrspos.size());
    ctx.restrspos = buffer_from_vector(parsed.restrspos, run_gpu);
    ctx.n_restrangs = static_cast<int>(parsed.restrangs.size());
    ctx.restrangs = buffer_from_vector(parsed.restrangs, run_gpu);
    ctx.n_restrdists = static_cast<int>(parsed.restrdists.size());
    ctx.restrdists = buffer_from_vector(parsed.restrdists, run_gpu);
    ctx.n_restrseqs = static_cast<int>(parsed.restrseqs.size());
    ctx.restrseqs = buffer_from_vector(parsed.restrseqs, run_gpu);
    ctx.n_restrwalls = static_cast<int>(parsed.restrwalls.size());
    ctx.restrwalls = buffer_from_vector(parsed.restrwalls, run_gpu);

    ctx.n_charges = static_cast<int>(parsed.charges.size());
    ctx.charges = buffer_from_vector(parsed.charges, run_gpu);
    ctx.n_ccharges = static_cast<int>(parsed.ccharges.size());
    ctx.ccharges = buffer_from_vector(parsed.ccharges, run_gpu);
    ctx.n_atypes = static_cast<int>(parsed.atypes.size());
    ctx.atypes = buffer_from_vector(parsed.atypes, run_gpu);
    ctx.n_catypes = static_cast<int>(parsed.catypes.size());
    ctx.catypes = buffer_from_vector(parsed.catypes, run_gpu);

    ctx.excluded = bool_buffer_from_vector(parsed.excluded, run_gpu);
    ctx.n_excluded = 0;
    for (bool excluded : parsed.excluded) {
        if (excluded) ctx.n_excluded++;
    }

    ctx.n_molecules = static_cast<int>(parsed.molecules.size());
    ctx.molecules = parsed.molecules;
    ctx.charge_group_config = parsed.charge_groups.empty() ? charge_group_config_t{} : parsed.charge_groups.front();

    ctx.n_qatoms = static_cast<int>(parsed.q_atoms.size());
    ctx.q_atoms = parsed.q_atoms;
    ctx.n_qangcouples = static_cast<int>(parsed.q_angcouples.size());
    ctx.q_angcouples = parsed.q_angcouples;
    ctx.n_qcatypes = static_cast<int>(parsed.q_catypes.size());
    ctx.q_catypes = parsed.q_catypes;
    ctx.q_atypes = parsed.q_atypes;
    ctx.q_charges = parsed.q_charges;

    ctx.n_qangles = per_lambda_count(parsed.q_angles.size(), ctx.n_lambdas);
    ctx.q_angles = parsed.q_angles;
    ctx.n_qcangles = static_cast<int>(parsed.q_cangles.size());
    ctx.q_cangles = parsed.q_cangles;
    ctx.n_qbonds = per_lambda_count(parsed.q_bonds.size(), ctx.n_lambdas);
    ctx.q_bonds = parsed.q_bonds;
    ctx.n_qcbonds = static_cast<int>(parsed.q_cbonds.size());
    ctx.q_cbonds = parsed.q_cbonds;
    ctx.n_qimpropers = per_lambda_count(parsed.q_impropers.size(), ctx.n_lambdas);
    ctx.q_impropers = parsed.q_impropers;
    ctx.n_qcimpropers = static_cast<int>(parsed.q_cimpropers.size());
    ctx.q_cimpropers = parsed.q_cimpropers;
    ctx.n_qtorsions = per_lambda_count(parsed.q_torsions.size(), ctx.n_lambdas);
    ctx.q_torsions = parsed.q_torsions;
    ctx.n_qctorsions = static_cast<int>(parsed.q_ctorsions.size());
    ctx.q_ctorsions = parsed.q_ctorsions;
    ctx.n_qoffdiags = static_cast<int>(parsed.q_offdiags.size());
    ctx.q_offdiags = parsed.q_offdiags;
    ctx.n_qimprcouples = static_cast<int>(parsed.q_imprcouples.size());
    ctx.q_imprcouples = parsed.q_imprcouples;
    ctx.n_qsoftpairs = static_cast<int>(parsed.q_softpairs.size());
    ctx.q_softpairs = parsed.q_softpairs;
    ctx.n_qtorcouples = static_cast<int>(parsed.q_torcouples.size());
    ctx.q_torcouples = parsed.q_torcouples;
    ctx.n_qelscales = per_lambda_count(parsed.q_elscales.size(), ctx.n_lambdas);
    ctx.q_elscales = buffer_from_vector(parsed.q_elscales, run_gpu);
    ctx.n_qexclpairs = per_lambda_count(parsed.q_exclpairs.size(), ctx.n_lambdas);
    ctx.q_exclpairs = parsed.q_exclpairs;
    ctx.n_qshakes = per_lambda_count(parsed.q_shakes.size(), ctx.n_lambdas);
    ctx.q_shakes = parsed.q_shakes;
    ctx.n_qsoftcores = per_lambda_count(parsed.q_softcores.size(), ctx.n_lambdas);
    ctx.q_softcores = parsed.q_softcores;

    set_lj_matrix(ctx, parsed);
}
}  // namespace

Context::Context() {
    if (current_ != nullptr) {
        fatal("Only one active Context instance is supported.");
        std::abort();
    }
    current_ = this;
}

Context::~Context() {
    if (current_ == this) {
        current_ = nullptr;
    }
}

Context& Context::instance() {
    if (current_ == nullptr) {
        fatal("Context::instance() called without an active Context.");
        std::abort();
    }
    return *current_;
}

void Context::cuda_reset_energies() {
    cudaMemset(dvelocities->gpu_data_p, 0, sizeof(dvel_t) * n_atoms);
    cudaMemset(EQ_restraint->gpu_data_p, 0, sizeof(E_restraint_t) * n_lambdas);
}

void Context::init_data_from_files() {
    std::unique_ptr<BaseParser> parser;
    if (command_info.input_mode == CommandInputMode::Csv) {
        parser = std::make_unique<CsvParser>(command_info.csv_dir);
    } else {
        parser = std::make_unique<InpParser>(command_info.input_file);
    }
    ParseResult parser_result = parser->parse();
    apply_parse_result(*this, parser_result);

    if (n_lambdas > 2) {
        fatal("More than 2 states not supported on GPU architecture. Exiting...");
    }
}

void Context::preprocess_data() {
    dt = time_unit * md.stepsize;
    tau_T = time_unit * md.bath_coupling;
    n_waters = (n_atoms - n_atoms_solute) / 3;
    init_inv_mass();
    exclude_qatom_definitions();

    // Shake constraints, need to be initialized before last part of shrink_topology
    if (md.charge_groups) {
        init_pshells_from_charge_groups();
    } else {
        init_pshells();
    }

    init_water_shell_parameters();
    init_shake();

    // Now remove shaken bonds
    exclude_shaken_definitions();

    preprocess_vdw_rule_parameters(*this);
    upload_preprocessed_topology(*this);

    init_unified_atom_parameters();

    // Init random seed from MD file
    srand(md.random_seed);

    init_patoms();
    init_velocities();

    dvelocities = std::make_unique<HostDeviceBuffer<dvel_t>>(n_atoms, true, command_info.requested_gpu);
    xcoords = std::make_unique<HostDeviceBuffer<coord_t>>(n_atoms, true, command_info.requested_gpu);

    if (n_waters > 0) {
        init_water_sphere();
        init_wshells();
    }

    // Init energy
    EQ_total.assign(n_lambdas, {});
    EQ_bond.assign(n_lambdas, {});
    EQ_nonbond_qq.assign(n_lambdas, {});
    EQ_nonbond_qp.assign(n_lambdas, {});
    EQ_nonbond_qw.assign(n_lambdas, {});
    EQ_nonbond_qx.assign(n_lambdas, {});
    EQ_restraint = std::make_unique<HostDeviceBuffer<E_restraint_t>>(n_lambdas, true, command_info.requested_gpu);

    if (md.has_initial_temperature && md.random_seed > 0) {
        initial_shaking();
        stop_cm_translation();
        if (command_info.requested_gpu) {
            coords->upload();
            velocities->upload();
            xcoords->upload();
        }
    }

    finalize_ngbrs14();

    init_atoms_list();
    initialize_charge_tables();
    initialize_catype_tables();
}

void Context::init() {
    init_data_from_files();
    preprocess_data();
}
