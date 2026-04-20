#include "cuda_context.cuh"
#include "context.h"
#include "common/include/constants.h"
#include "common/include/vdw_rules.h"

void CudaContext::init() {
    auto& host = Context::instance();

    check_cudaMalloc((void**)&d_EQ_nonbond_qq, sizeof(E_nonbonded_t) * host.n_lambdas);
    check_cudaMalloc((void**)&d_lambdas, sizeof(double) * host.n_lambdas);

    check_cudaMalloc((void**)&d_wshells, host.n_shells * sizeof(shell_t));
    check_cudaMalloc((void**)&d_shell, sizeof(bool) * host.n_atoms);

    check_cudaMalloc((void**)&d_EQ_restraint, sizeof(E_restraint_t) * host.n_lambdas);

    check_cudaMalloc((void**)&d_unified_ccharges, sizeof(ccharge_t) * host.unified_ccharges.size());
    check_cudaMalloc((void**)&d_unified_catypes, sizeof(catype_t) * host.unified_catypes.size());

    initialize_atom_lists_host();
    initialize_ngbrs14_host();
    initialize_charge_tables_host();
    initialize_catype_tables_host();

    sync_all_to_device();
}

void CudaContext::sync_all_to_device() {
    auto& host = Context::instance();

    host.coords->upload();
    host.dvelocities->upload();
    host.velocities->upload();
    host.LJ_matrix->upload();
    host.xcoords->upload();
    host.winv->upload();
    host.mol_n_shakes->upload();
    host.shake_bonds->upload();
    host.heavy->upload();
    host.excluded->upload();
    host.q_elscales->upload();
    sync_array_to_device<E_nonbonded_t>(d_EQ_nonbond_qq, host.EQ_nonbond_qq.data(), host.n_lambdas);
    sync_array_to_device<double>(d_lambdas, host.lambdas.data(), host.n_lambdas);
    sync_array_to_device<shell_t>(d_wshells, host.wshells.data(), host.n_shells);

    sync_array_to_device<bool>(d_shell, host.shell.get(), host.n_atoms);

    sync_array_to_device<E_restraint_t>(d_EQ_restraint, host.EQ_restraint.data(), host.n_lambdas);

    sync_array_to_device<ccharge_t>(d_unified_ccharges, host.unified_ccharges.data(), host.unified_ccharges.size());
    sync_array_to_device<catype_t>(d_unified_catypes, host.unified_catypes.data(), host.unified_catypes.size());
}

void CudaContext::sync_all_to_host() {
    auto& host = Context::instance();

    host.coords->download();
    host.dvelocities->download();
    host.velocities->download();
    host.LJ_matrix->download();
    host.xcoords->download();
    host.winv->download();
    host.mol_n_shakes->download();
    host.shake_bonds->download();
    host.excluded->download();
    host.q_elscales->download();
    sync_array_to_host<E_nonbonded_t>(host.EQ_nonbond_qq.data(), d_EQ_nonbond_qq, host.n_lambdas);
    sync_array_to_host<double>(host.lambdas.data(), d_lambdas, host.n_lambdas);
    sync_array_to_host<shell_t>(host.wshells.data(), d_wshells, host.n_shells);
    sync_array_to_host<bool>(host.shell.get(), d_shell, host.n_atoms);

    sync_array_to_host<E_restraint_t>(host.EQ_restraint.data(), d_EQ_restraint, host.n_lambdas);

    sync_array_to_host<ccharge_t>(host.unified_ccharges.data(), d_unified_ccharges, host.unified_ccharges.size());
    sync_array_to_host<catype_t>(host.unified_catypes.data(), d_unified_catypes, host.unified_catypes.size());
}

void CudaContext::free() {
    cudaFree(d_EQ_nonbond_qq);
    cudaFree(d_lambdas);

    cudaFree(d_wshells);
    cudaFree(d_shell);

    cudaFree(d_EQ_restraint);

    cudaFree(d_charge_table_all);
    cudaFree(d_charge_pair_products);
    cudaFree(d_catype_table_all);
    cudaFree(d_catype_pair_params);
    cudaFree(d_unified_ccharges);
    cudaFree(d_unified_catypes);
    cudaFree(d_ngbrs_14);
    cudaFree(d_p_charge_types);
    cudaFree(d_w_charge_types);
    cudaFree(d_q_charge_types);
    cudaFree(d_p_catype_types);
    cudaFree(d_w_catype_types);
    cudaFree(d_q_catype_types);
    cudaFree(d_p_atoms_list);
    cudaFree(d_w_atoms_list);
    cudaFree(d_q_atoms_list);

    d_charge_table_all = nullptr;
    d_charge_pair_products = nullptr;
    d_catype_table_all = nullptr;
    d_catype_pair_params = nullptr;
    d_unified_ccharges = nullptr;
    d_unified_catypes = nullptr;
    d_ngbrs_14 = nullptr;
    n_charge_types = 0;
    zero_charge_type = -1;
    n_catype_types = 0;
    zero_catype_type = -1;
    d_p_charge_types = nullptr;
    d_w_charge_types = nullptr;
    d_q_charge_types = nullptr;
    d_p_catype_types = nullptr;
    d_w_catype_types = nullptr;
    d_q_catype_types = nullptr;
    d_p_atoms_list = nullptr;
    d_w_atoms_list = nullptr;
    d_q_atoms_list = nullptr;
}

void CudaContext::reset_energies() {
    auto& host = Context::instance();
    cudaMemset(host.dvelocities->gpu_data_p, 0, sizeof(dvel_t) * host.n_atoms);
    cudaMemset(d_EQ_nonbond_qq, 0, sizeof(E_nonbonded_t) * host.n_lambdas);
    cudaMemset(d_EQ_restraint, 0, sizeof(E_restraint_t) * host.n_lambdas);
}

void CudaContext::initialize_atom_lists_host() {
    auto& host = Context::instance();
    auto *excluded = host.excluded->cpu_data_p;

    h_p_atoms_list.clear();
    for (int i = 0; i < host.n_patoms; i++) {
        int id = host.p_atoms[i];
        if (!excluded[id]) {
            h_p_atoms_list.push_back(id);
        }
    }

    h_w_atoms_list.clear();
    for (int i = host.n_atoms_solute; i < host.n_atoms; i++) {
        int id = i;
        if (!excluded[id]) {
            h_w_atoms_list.push_back(id);
        }
    }
    printf("Number of water atoms: %d, number of all water atoms: %d\n", (int)h_w_atoms_list.size(), host.n_atoms - host.n_atoms_solute);

    h_q_atoms_list.clear();
    for (int i = 0; i < host.n_qatoms; i++) {
        int id = host.q_atoms[i];
        if (!excluded[id]) {
            h_q_atoms_list.push_back(id);
        }
    }

    check_cudaMalloc((void**)&d_p_atoms_list, sizeof(int) * h_p_atoms_list.size());
    cudaMemcpy(d_p_atoms_list, h_p_atoms_list.data(), sizeof(int) * h_p_atoms_list.size(), cudaMemcpyHostToDevice);
    check_cudaMalloc((void**)&d_w_atoms_list, sizeof(int) * h_w_atoms_list.size());
    cudaMemcpy(d_w_atoms_list, h_w_atoms_list.data(), sizeof(int) * h_w_atoms_list.size(), cudaMemcpyHostToDevice);
    check_cudaMalloc((void**)&d_q_atoms_list, sizeof(int) * h_q_atoms_list.size());
    cudaMemcpy(d_q_atoms_list, h_q_atoms_list.data(), sizeof(int) * h_q_atoms_list.size(), cudaMemcpyHostToDevice);
}

void CudaContext::initialize_ngbrs14_host() {
    auto& host = Context::instance();

    if (!host.ngbrs_14.empty()) {
        check_cudaMalloc((void**)&d_ngbrs_14, sizeof(int3) * host.ngbrs_14.size());
        cudaMemcpy(d_ngbrs_14, host.ngbrs_14.data(), sizeof(int3) * host.ngbrs_14.size(), cudaMemcpyHostToDevice);
    }
}

void CudaContext::initialize_charge_tables_host() {
    auto& host = Context::instance();
    auto *excluded = host.excluded->cpu_data_p;
    auto &charges = host.charges->cpu_data_p;
    auto &ccharges = host.ccharges->cpu_data_p;

    std::map<double, int> charge_to_type_host;
    std::vector<ccharge_t> h_charge_table_all;  // h_charge_table_all[charge type] = (code, charge)

    auto add_charge = [&](double charge) -> int {
        if (charge_to_type_host.count(charge) == 0) {
            int sz = h_charge_table_all.size();
            ccharge_t new_ccharge;
            new_ccharge.code = sz;
            new_ccharge.charge = charge;
            h_charge_table_all.push_back(new_ccharge);
            charge_to_type_host[charge] = sz;
        }
        return charge_to_type_host[charge];
    };

    for (int i = 0; i < host.n_ccharges; i++) {
        double charge = ccharges[i].charge;
        add_charge(charge);
    }
    for (int state = 0; state < host.n_lambdas; state++) {
        for (int i = 0; i < host.n_qatoms; i++) {
            double charge = host.q_charges[i + host.n_qatoms * state].charge;
            add_charge(charge);
            charge *= host.lambdas[state];
            add_charge(charge);
        }
    }

    zero_charge_type = add_charge(0.0);
    n_charge_types = static_cast<int>(h_charge_table_all.size());

    std::vector<double> h_charge_pair_products(n_charge_types * n_charge_types);
    for (int i = 0; i < n_charge_types; i++) {
        for (int j = 0; j < n_charge_types; j++) {
            h_charge_pair_products[i * n_charge_types + j] = h_charge_table_all[i].charge * h_charge_table_all[j].charge;
        }
    }

    std::vector<int> p_charge_types(h_p_atoms_list.size());
    for (int i = 0; i < static_cast<int>(h_p_atoms_list.size()); i++) {
        int id = h_p_atoms_list[i];
        double charge = ccharges[charges[id].code - 1].charge;
        p_charge_types[i] = charge_to_type_host[charge];
    }

    int h_q_atoms_size = h_q_atoms_list.size();
    std::vector<int> q_charge_types(h_q_atoms_size * host.n_lambdas);
    std::map<int, int> q_idx;
    for (int i = 0; i < host.n_qatoms; i++) {
        int id = host.q_atoms[i];
        if (!excluded[id]) {
            q_idx[id] = i;
        }
    }

    for (int state = 0; state < host.n_lambdas; state++) {
        for (int i = 0; i < h_q_atoms_size; i++) {
            int id = h_q_atoms_list[i];
            double charge = host.q_charges[q_idx[id] + host.n_qatoms * state].charge;
            q_charge_types[state * h_q_atoms_size + i] = charge_to_type_host[charge];
        }
    }

    int h_w_atoms_size = h_w_atoms_list.size();
    std::vector<int> w_charge_types(h_w_atoms_size);
    for (int i = 0; i < h_w_atoms_size; i++) {
        int id = h_w_atoms_list[i];
        double charge = ccharges[charges[id].code - 1].charge;
        w_charge_types[i] = charge_to_type_host[charge];
    }

    check_cudaMalloc((void**)&d_charge_table_all, sizeof(ccharge_t) * h_charge_table_all.size());
    cudaMemcpy(d_charge_table_all, h_charge_table_all.data(), sizeof(ccharge_t) * h_charge_table_all.size(), cudaMemcpyHostToDevice);
    check_cudaMalloc((void**)&d_charge_pair_products, sizeof(double) * h_charge_pair_products.size());
    cudaMemcpy(d_charge_pair_products, h_charge_pair_products.data(), sizeof(double) * h_charge_pair_products.size(), cudaMemcpyHostToDevice);
    check_cudaMalloc((void**)&d_p_charge_types, sizeof(int) * p_charge_types.size());
    cudaMemcpy(d_p_charge_types, p_charge_types.data(), sizeof(int) * p_charge_types.size(), cudaMemcpyHostToDevice);
    check_cudaMalloc((void**)&d_q_charge_types, sizeof(int) * q_charge_types.size());
    cudaMemcpy(d_q_charge_types, q_charge_types.data(), sizeof(int) * q_charge_types.size(), cudaMemcpyHostToDevice);
    check_cudaMalloc((void**)&d_w_charge_types, sizeof(int) * w_charge_types.size());
    cudaMemcpy(d_w_charge_types, w_charge_types.data(), sizeof(int) * w_charge_types.size(), cudaMemcpyHostToDevice);
}

void CudaContext::initialize_catype_tables_host() {
    auto& host = Context::instance();
    auto *excluded = host.excluded->cpu_data_p;
    auto &atypes = host.atypes->cpu_data_p;
    auto &catypes = host.catypes->cpu_data_p;

    std::vector<catype_t> h_catype_table_all;  // h_catype_table_all[catype code] = catype_t
    catype_to_type_host.clear();

    auto add_catype = [&](catype_t catype) -> int {
        auto key = get_catype_key(catype);
        if (catype_to_type_host.count(key) == 0) {
            int sz = h_catype_table_all.size();
            catype.code = sz;
            h_catype_table_all.push_back(catype);
            catype_to_type_host[key] = sz;
        }
        return catype_to_type_host[key];
    };

    for (int i = 0; i < host.n_catypes; i++) {
        add_catype(catypes[i]);
    }

    for (int state = 0; state < host.n_lambdas; state++) {
        for (int i = 0; i < host.n_qatoms; i++) {
            const atype_t& qat = host.q_atypes[i + host.n_qatoms * state];
            add_catype(host.q_catypes[qat.code - 1]);
        }
    }

    catype_t zero_catype = {};
    zero_catype_type = add_catype(zero_catype);
    n_catype_types = static_cast<int>(h_catype_table_all.size());

    std::vector<vdw_pair_param_t> h_catype_pair_params(n_catype_types * n_catype_types);
    for (int i = 0; i < n_catype_types; i++) {
        for (int j = 0; j < n_catype_types; j++) {
            const catype_t& ci = h_catype_table_all[i];
            const catype_t& cj = h_catype_table_all[j];
            vdw_pair_param_t pair_param = {};
            if (host.topo.vdw_rule == VDW_GEOMETRIC) {
                calc_vdw_geometric(ci.aii_normal, cj.aii_normal, ci.bii_normal, cj.bii_normal, 1.0, &pair_param.a, &pair_param.b);
            } else {
                calc_vdw_arithmetic(ci.aii_normal, cj.aii_normal, ci.bii_normal, cj.bii_normal, 1.0, &pair_param.a, &pair_param.b);
            }
            h_catype_pair_params[i * n_catype_types + j] = pair_param;
        }
    }

    std::vector<int> p_catype_types(h_p_atoms_list.size());
    for (int i = 0; i < static_cast<int>(h_p_atoms_list.size()); i++) {
        int id = h_p_atoms_list[i];
        auto catype = catypes[atypes[id].code - 1];
        auto key = get_catype_key(catype);
        p_catype_types[i] = catype_to_type_host[key];
    }

    int h_q_atoms_size = h_q_atoms_list.size();
    std::vector<int> q_catype_types(h_q_atoms_size * host.n_lambdas);
    std::map<int, int> q_idx;
    for (int i = 0; i < host.n_qatoms; i++) {
        int id = host.q_atoms[i];
        if (!excluded[id]) {
            q_idx[id] = i;
        }
    }


    for (int state = 0; state < host.n_lambdas; state++) {
        for (int i = 0; i < h_q_atoms_size; i++) {
            int id = h_q_atoms_list[i];
            const atype_t& qat = host.q_atypes[q_idx[id] + host.n_qatoms * state];
            auto key_normal = get_catype_key(host.q_catypes[qat.code - 1]);
            q_catype_types[state * h_q_atoms_size + i] = catype_to_type_host[key_normal];
        }
    }

    int h_w_atoms_size = h_w_atoms_list.size();
    std::vector<int> w_catype_types(h_w_atoms_size);
    for (int i = 0; i < h_w_atoms_size; i++) {
        int id = h_w_atoms_list[i];
        auto catype = catypes[atypes[id].code - 1];
        auto key = get_catype_key(catype);
        w_catype_types[i] = catype_to_type_host[key];
        // printf("Water atom %d: catype code %d mapped to %d, data: m=%f, aii_normal=%f, bii_normal=%f, aii_polar=%f, bii_polar=%f, aii_1_4=%f, bii_1_4=%f\n", id, atypes[id].code, w_catype_types[i], catype.m, catype.aii_normal, catype.bii_normal, catype.aii_polar, catype.bii_polar, catype.aii_1_4, catype.bii_1_4);
    }
    printf("Total water atom number: %d, w_catype_types size: %lu\n", h_w_atoms_size, w_catype_types.size());

    check_cudaMalloc((void**)&d_catype_table_all, sizeof(catype_t) * h_catype_table_all.size());
    cudaMemcpy(d_catype_table_all, h_catype_table_all.data(), sizeof(catype_t) * h_catype_table_all.size(), cudaMemcpyHostToDevice);
    check_cudaMalloc((void**)&d_catype_pair_params, sizeof(vdw_pair_param_t) * h_catype_pair_params.size());
    cudaMemcpy(d_catype_pair_params, h_catype_pair_params.data(), sizeof(vdw_pair_param_t) * h_catype_pair_params.size(), cudaMemcpyHostToDevice);
    check_cudaMalloc((void**)&d_p_catype_types, sizeof(int) * p_catype_types.size());
    cudaMemcpy(d_p_catype_types, p_catype_types.data(), sizeof(int) * p_catype_types.size(), cudaMemcpyHostToDevice);
    check_cudaMalloc((void**)&d_q_catype_types, sizeof(int) * q_catype_types.size());
    cudaMemcpy(d_q_catype_types, q_catype_types.data(), sizeof(int) * q_catype_types.size(), cudaMemcpyHostToDevice);
    check_cudaMalloc((void**)&d_w_catype_types, sizeof(int) * w_catype_types.size());
    cudaMemcpy(d_w_catype_types, w_catype_types.data(), sizeof(int) * w_catype_types.size(), cudaMemcpyHostToDevice);
}
