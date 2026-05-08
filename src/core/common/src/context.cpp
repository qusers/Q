#include "context.h"

#include <cstdio>
#include <cstdlib>

#include "constants.h"
#include "init.h"
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
        std::fprintf(stderr, ">>> FATAL: Only one active Context instance is supported.\n");
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
        std::fprintf(stderr, ">>> FATAL: Context::instance() called without an active Context.\n");
        std::abort();
    }
    return *current_;
}

void Context::cuda_reset_energies() {
    cudaMemset(dvelocities->gpu_data_p, 0, sizeof(dvel_t) * n_atoms);
    cudaMemset(EQ_restraint->gpu_data_p, 0, sizeof(E_restraint_t) * n_lambdas);
}

void Context::init_data_from_files() {
    if (q_input_mode == Q_INPUT_NATIVE) {
        if (n_lambdas > 2) {
            printf(">>> FATAL: More than 2 states not supported on GPU architecture. Exiting...\n");
            exit(EXIT_FAILURE);
        }
        return;
    }

    parse_md("md.csv");
    if (n_lambdas > 2) {
        printf(">>> FATAL: More than 2 states not supported on GPU architecture. Exiting...\n");
        exit(EXIT_FAILURE);
    }
    parse_topo("topo.csv");
    parse_angles("angles.csv");
    parse_atypes("atypes.csv");
    parse_bonds("bonds.csv");
    parse_cangles("cangles.csv");
    parse_catypes("catypes.csv");
    parse_cbonds("cbonds.csv");
    parse_ccharges("ccharges.csv");
    parse_charges("charges.csv");
    parse_cimpropers("cimpropers.csv");
    parse_coords("coords.csv");
    parse_icoords("i_coords.csv");
    parse_ctorsions("ctorsions.csv");
    parse_excluded("excluded.csv");
    parse_molecules("molecules.csv");
    parse_impropers("impropers.csv");
    parse_torsions("torsions.csv");
    parse_LJ_matrix();
    parse_qatoms("q_atoms.csv");
    parse_qangcouples("q_angcouples.csv");
    parse_qcangles("q_cangles.csv");
    parse_qcatypes("q_catypes.csv");
    parse_qcbonds("q_cbonds.csv");
    parse_qcimpropers("q_cimpropers.csv");
    parse_qctorsions("q_ctorsions.csv");
    parse_qoffdiags("q_offdiags.csv");
    parse_qimprcouples("q_imprcouples.csv");
    parse_qsoftpairs("q_softpairs.csv");
    parse_qtorcouples("q_torcouples.csv");
    parse_qangles("q_angles.csv");
    parse_qatypes("q_atypes.csv");
    parse_qbonds("q_bonds.csv");
    parse_qcharges("q_charges.csv");
    parse_qelscales("q_elscales.csv");
    parse_qexclpairs("q_exclpairs.csv");
    parse_qimpropers("q_impropers.csv");
    parse_qshakes("q_shakes.csv");
    parse_qsoftcores("q_softcores.csv");
    parse_qtorsions("q_torsions.csv");
    if (md.charge_groups) {
        parse_charge_groups("charge_groups.csv");
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

    dvelocities = std::make_unique<HostDeviceBuffer<dvel_t>>(n_atoms, true, run_gpu);
    xcoords = std::make_unique<HostDeviceBuffer<coord_t>>(n_atoms, true, run_gpu);

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
    EQ_restraint = std::make_unique<HostDeviceBuffer<E_restraint_t>>(n_lambdas, true, run_gpu);

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
    if (run_gpu) {
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
