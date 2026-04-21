#include <algorithm>
#include <map>
#include <vector>

#include "context.h"
#include "common/include/constants.h"
#include "common/include/vdw_rules.h"

namespace {
template <typename T>
std::unique_ptr<HostDeviceBuffer<T>> make_host_device_buffer_from_vector(
    const std::vector<T>& src,
    bool run_gpu) {
    auto buffer = std::make_unique<HostDeviceBuffer<T>>(src.size(), true, run_gpu);
    if (!src.empty()) {
        std::copy(src.begin(), src.end(), buffer->cpu_data_p);
    }
    return buffer;
}
}  // namespace

void Context::cuda_initialize_helpers() {
    cuda_initialize_atom_lists_host();
    cuda_initialize_ngbrs14_host();
    cuda_initialize_charge_tables_host();
    cuda_initialize_catype_tables_host();
    cuda_sync_all_to_device();
}

void Context::cuda_sync_all_to_device() {
    coords->upload();
    dvelocities->upload();
    velocities->upload();
    LJ_matrix->upload();
    xcoords->upload();
    winv->upload();
    mol_n_shakes->upload();
    shake_bonds->upload();
    heavy->upload();
    excluded->upload();
    q_elscales->upload();
    shell->upload();
    wshells->upload();
    lambdas->upload();
    EQ_restraint->upload();
    unified_ccharges->upload();
    unified_catypes->upload();
    ngbrs_14->upload();
    p_atoms_list->upload();
    w_atoms_list->upload();
    q_atoms_list->upload();
    charge_table_all->upload();
    charge_pair_products->upload();
    p_charge_types->upload();
    w_charge_types->upload();
    q_charge_types->upload();
    catype_table_all->upload();
    catype_pair_params->upload();
    p_catype_types->upload();
    w_catype_types->upload();
    q_catype_types->upload();
}

void Context::cuda_sync_all_to_host() {
    coords->download();
    dvelocities->download();
    velocities->download();
    LJ_matrix->download();
    xcoords->download();
    winv->download();
    mol_n_shakes->download();
    shake_bonds->download();
    excluded->download();
    q_elscales->download();
    shell->download();
    wshells->download();
    lambdas->download();
    EQ_restraint->download();
    unified_ccharges->download();
    unified_catypes->download();
    ngbrs_14->download();
    p_atoms_list->download();
    w_atoms_list->download();
    q_atoms_list->download();
    charge_table_all->download();
    charge_pair_products->download();
    p_charge_types->download();
    w_charge_types->download();
    q_charge_types->download();
    catype_table_all->download();
    catype_pair_params->download();
    p_catype_types->download();
    w_catype_types->download();
    q_catype_types->download();
}

void Context::cuda_free_helpers() {
    unified_ccharges.reset();
    unified_catypes.reset();
    ngbrs_14.reset();
    p_atoms_list.reset();
    w_atoms_list.reset();
    q_atoms_list.reset();
    charge_table_all.reset();
    charge_pair_products.reset();
    p_charge_types.reset();
    w_charge_types.reset();
    q_charge_types.reset();
    catype_table_all.reset();
    catype_pair_params.reset();
    p_catype_types.reset();
    w_catype_types.reset();
    q_catype_types.reset();
    catype_to_type_host.clear();
    n_charge_types = 0;
    zero_charge_type = -1;
    n_catype_types = 0;
    zero_catype_type = -1;
}

void Context::cuda_reset_energies() {
    cudaMemset(dvelocities->gpu_data_p, 0, sizeof(dvel_t) * n_atoms);
    cudaMemset(EQ_restraint->gpu_data_p, 0, sizeof(E_restraint_t) * n_lambdas);
}

void Context::cuda_initialize_atom_lists_host() {
    auto *excluded = this->excluded->cpu_data_p;

    std::vector<int> h_p_atoms_list;
    h_p_atoms_list.reserve(n_patoms);
    for (int i = 0; i < n_patoms; i++) {
        int id = p_atoms[i];
        if (!excluded[id]) {
            h_p_atoms_list.push_back(id);
        }
    }

    std::vector<int> h_w_atoms_list;
    h_w_atoms_list.reserve(n_atoms - n_atoms_solute);
    for (int i = n_atoms_solute; i < n_atoms; i++) {
        if (!excluded[i]) {
            h_w_atoms_list.push_back(i);
        }
    }
    printf("Number of water atoms: %d, number of all water atoms: %d\n", (int)h_w_atoms_list.size(), n_atoms - n_atoms_solute);

    std::vector<int> h_q_atoms_list;
    h_q_atoms_list.reserve(n_qatoms);
    for (int i = 0; i < n_qatoms; i++) {
        int id = q_atoms[i];
        if (!excluded[id]) {
            h_q_atoms_list.push_back(id);
        }
    }

    p_atoms_list = make_host_device_buffer_from_vector(h_p_atoms_list, run_gpu);
    w_atoms_list = make_host_device_buffer_from_vector(h_w_atoms_list, run_gpu);
    q_atoms_list = make_host_device_buffer_from_vector(h_q_atoms_list, run_gpu);
}

void Context::cuda_initialize_ngbrs14_host() {
    if (ngbrs_14 == nullptr) {
        ngbrs_14 = std::make_unique<HostDeviceBuffer<int3>>(0, true, run_gpu);
    }
}

void Context::cuda_initialize_charge_tables_host() {
    auto *excluded = this->excluded->cpu_data_p;
    auto &charges = this->charges->cpu_data_p;
    auto &ccharges = this->ccharges->cpu_data_p;
    auto *lambda_values = lambdas->cpu_data_p;
    auto *p_atoms_cpu = p_atoms_list->cpu_data_p;
    auto *w_atoms_cpu = w_atoms_list->cpu_data_p;
    auto *q_atoms_cpu = q_atoms_list->cpu_data_p;

    std::map<double, int> charge_to_type_host;
    std::vector<ccharge_t> h_charge_table_all;

    auto add_charge = [&](double charge) -> int {
        if (charge_to_type_host.count(charge) == 0) {
            int sz = static_cast<int>(h_charge_table_all.size());
            ccharge_t new_ccharge = {};
            new_ccharge.code = sz;
            new_ccharge.charge = charge;
            h_charge_table_all.push_back(new_ccharge);
            charge_to_type_host[charge] = sz;
        }
        return charge_to_type_host[charge];
    };

    for (int i = 0; i < n_ccharges; i++) {
        add_charge(ccharges[i].charge);
    }
    for (int state = 0; state < n_lambdas; state++) {
        for (int i = 0; i < n_qatoms; i++) {
            double charge = q_charges[i + n_qatoms * state].charge;
            add_charge(charge);
            add_charge(charge * lambda_values[state]);
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

    std::vector<int> p_charge_types_cpu(p_atoms_list->length);
    for (int i = 0; i < static_cast<int>(p_atoms_list->length); i++) {
        const int id = p_atoms_cpu[i];
        const double charge = ccharges[charges[id].code - 1].charge;
        p_charge_types_cpu[i] = charge_to_type_host[charge];
    }

    std::vector<int> q_charge_types_cpu(q_atoms_list->length * n_lambdas);
    std::map<int, int> q_idx;
    for (int i = 0; i < n_qatoms; i++) {
        int id = q_atoms[i];
        if (!excluded[id]) {
            q_idx[id] = i;
        }
    }

    for (int state = 0; state < n_lambdas; state++) {
        for (int i = 0; i < static_cast<int>(q_atoms_list->length); i++) {
            const int id = q_atoms_cpu[i];
            const double charge = q_charges[q_idx[id] + n_qatoms * state].charge;
            q_charge_types_cpu[state * q_atoms_list->length + i] = charge_to_type_host[charge];
        }
    }

    std::vector<int> w_charge_types_cpu(w_atoms_list->length);
    for (int i = 0; i < static_cast<int>(w_atoms_list->length); i++) {
        const int id = w_atoms_cpu[i];
        const double charge = ccharges[charges[id].code - 1].charge;
        w_charge_types_cpu[i] = charge_to_type_host[charge];
    }

    charge_table_all = make_host_device_buffer_from_vector(h_charge_table_all, run_gpu);
    charge_pair_products = make_host_device_buffer_from_vector(h_charge_pair_products, run_gpu);
    p_charge_types = make_host_device_buffer_from_vector(p_charge_types_cpu, run_gpu);
    q_charge_types = make_host_device_buffer_from_vector(q_charge_types_cpu, run_gpu);
    w_charge_types = make_host_device_buffer_from_vector(w_charge_types_cpu, run_gpu);
}

void Context::cuda_initialize_catype_tables_host() {
    auto *excluded = this->excluded->cpu_data_p;
    auto &atypes = this->atypes->cpu_data_p;
    auto &catypes = this->catypes->cpu_data_p;
    auto *p_atoms_cpu = p_atoms_list->cpu_data_p;
    auto *w_atoms_cpu = w_atoms_list->cpu_data_p;
    auto *q_atoms_cpu = q_atoms_list->cpu_data_p;

    std::vector<catype_t> h_catype_table_all;
    catype_to_type_host.clear();

    auto add_catype = [&](catype_t catype) -> int {
        const std::array<double, 4> key = {
            catype.aii_normal,
            catype.bii_normal,
            catype.aii_1_4,
            catype.bii_1_4};
        if (catype_to_type_host.count(key) == 0) {
            int sz = static_cast<int>(h_catype_table_all.size());
            catype.code = sz;
            h_catype_table_all.push_back(catype);
            catype_to_type_host[key] = sz;
        }
        return catype_to_type_host[key];
    };

    for (int i = 0; i < n_catypes; i++) {
        add_catype(catypes[i]);
    }

    for (int state = 0; state < n_lambdas; state++) {
        for (int i = 0; i < n_qatoms; i++) {
            const atype_t& qat = q_atypes[i + n_qatoms * state];
            add_catype(q_catypes[qat.code - 1]);
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
            if (topo.vdw_rule == VDW_GEOMETRIC) {
                calc_vdw_geometric(ci.aii_normal, cj.aii_normal, ci.bii_normal, cj.bii_normal, 1.0, &pair_param.a, &pair_param.b);
            } else {
                calc_vdw_arithmetic(ci.aii_normal, cj.aii_normal, ci.bii_normal, cj.bii_normal, 1.0, &pair_param.a, &pair_param.b);
            }
            h_catype_pair_params[i * n_catype_types + j] = pair_param;
        }
    }

    std::vector<int> p_catype_types_cpu(p_atoms_list->length);
    for (int i = 0; i < static_cast<int>(p_atoms_list->length); i++) {
        const int id = p_atoms_cpu[i];
        const catype_t catype = catypes[atypes[id].code - 1];
        const std::array<double, 4> key = {catype.aii_normal, catype.bii_normal, catype.aii_1_4, catype.bii_1_4};
        p_catype_types_cpu[i] = catype_to_type_host[key];
    }

    std::vector<int> q_catype_types_cpu(q_atoms_list->length * n_lambdas);
    std::map<int, int> q_idx;
    for (int i = 0; i < n_qatoms; i++) {
        int id = q_atoms[i];
        if (!excluded[id]) {
            q_idx[id] = i;
        }
    }

    for (int state = 0; state < n_lambdas; state++) {
        for (int i = 0; i < static_cast<int>(q_atoms_list->length); i++) {
            const int id = q_atoms_cpu[i];
            const atype_t& qat = q_atypes[q_idx[id] + n_qatoms * state];
            const catype_t& qcatype = q_catypes[qat.code - 1];
            const std::array<double, 4> key = {qcatype.aii_normal, qcatype.bii_normal, qcatype.aii_1_4, qcatype.bii_1_4};
            q_catype_types_cpu[state * q_atoms_list->length + i] = catype_to_type_host[key];
        }
    }

    std::vector<int> w_catype_types_cpu(w_atoms_list->length);
    for (int i = 0; i < static_cast<int>(w_atoms_list->length); i++) {
        const int id = w_atoms_cpu[i];
        const catype_t catype = catypes[atypes[id].code - 1];
        const std::array<double, 4> key = {catype.aii_normal, catype.bii_normal, catype.aii_1_4, catype.bii_1_4};
        w_catype_types_cpu[i] = catype_to_type_host[key];
    }
    printf("Total water atom number: %lu, w_catype_types size: %lu\n", w_atoms_list->length, w_catype_types_cpu.size());

    catype_table_all = make_host_device_buffer_from_vector(h_catype_table_all, run_gpu);
    catype_pair_params = make_host_device_buffer_from_vector(h_catype_pair_params, run_gpu);
    p_catype_types = make_host_device_buffer_from_vector(p_catype_types_cpu, run_gpu);
    q_catype_types = make_host_device_buffer_from_vector(q_catype_types_cpu, run_gpu);
    w_catype_types = make_host_device_buffer_from_vector(w_catype_types_cpu, run_gpu);
}
