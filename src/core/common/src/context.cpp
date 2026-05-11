#include "context.h"

#include <cstdio>
#include <cstdlib>

#include "constants.h"
#include "csv_parser.h"
#include "helpers.h"
#include "init.h"
#include "inp_parser.h"
#include "parse.h"
#include "q_input.h"

Context* Context::current_ = nullptr;

namespace {
template <typename T>
void upload_if_present(const std::unique_ptr<HostDeviceBuffer<T>>& buffer) {
    if (buffer) {
        buffer->upload();
    }
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
    auto parser_result = parser->parse();

    if (n_lambdas > 2) {
        fatal("More than 2 states not supported on GPU architecture. Exiting...");
        exit(EXIT_FAILURE);
    }
}

void Context::preprocess_data() {
    dt = time_unit * md.stepsize;
    tau_T = time_unit * md.bath_coupling;
    init_inv_mass();
    exclude_qatom_definitions();

    // Shake constraints, need to be initialized before last part of shrink_topology
    if (md.charge_groups) {
        init_pshells_from_charge_groups();
    } else {
        init_pshells();
    }

    init_shake();

    // Now remove shaken bonds
    exclude_shaken_definitions();

    init_unified_atom_parameters();

    // Init random seed from MD file
    srand(md.random_seed);

    init_patoms();
    init_velocities();

    dvelocities = std::make_unique<HostDeviceBuffer<dvel_t>>(n_atoms, true, command_info.requested_gpu);
    xcoords = std::make_unique<HostDeviceBuffer<coord_t>>(n_atoms, true, command_info.requested_gpu);

    n_waters = (n_atoms - n_atoms_solute) / 3;
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

    if (n_shake_constraints > 0 || (q_input_mode == Q_INPUT_NATIVE && q_native_fresh_start && md.random_seed > 0)) {
        initial_shaking();
        stop_cm_translation();
    }

    finalize_ngbrs14();

    init_atoms_list();
    initialize_charge_tables();
    initialize_catype_tables();
}

void Context::gpu_data_upload() {
    if (command_info.requested_gpu) {
        upload_if_present(coords_init);
        upload_if_present(coords);
        upload_if_present(velocities);
        upload_if_present(dvelocities);
        upload_if_present(xcoords);

        upload_if_present(angles);
        upload_if_present(cangles);
        upload_if_present(bonds);
        upload_if_present(cbonds);
        upload_if_present(torsions);
        upload_if_present(ctorsions);
        upload_if_present(impropers);
        upload_if_present(cimpropers);

        upload_if_present(restrspos);
        upload_if_present(restrangs);
        upload_if_present(restrdists);
        upload_if_present(restrseqs);
        upload_if_present(restrwalls);

        upload_if_present(charges);
        upload_if_present(ccharges);
        upload_if_present(atypes);
        upload_if_present(catypes);
        upload_if_present(unified_ccharges);
        upload_if_present(unified_catypes);
        upload_if_present(excluded);
        upload_if_present(heavy);
        upload_if_present(shell);
        upload_if_present(winv);

        upload_if_present(lambdas);
        upload_if_present(q_elscales);

        upload_if_present(LJ_matrix);
        upload_if_present(ngbrs_14);

        upload_if_present(mol_n_shakes);
        upload_if_present(shake_bonds);
        upload_if_present(wshells);

        upload_if_present(EQ_restraint);

        upload_if_present(p_atoms_list);
        upload_if_present(w_atoms_list);
        upload_if_present(q_atoms_list);
        upload_if_present(charge_pair_products);
        upload_if_present(p_charge_types);
        upload_if_present(w_charge_types);
        upload_if_present(q_charge_types);
        upload_if_present(catype_pair_params);
        upload_if_present(p_catype_types);
        upload_if_present(w_catype_types);
        upload_if_present(q_catype_types);
    }
}

void Context::init() {
    init_data_from_files();
    preprocess_data();
    gpu_data_upload();
}
