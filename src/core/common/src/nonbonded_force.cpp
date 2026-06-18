#include "nonbonded_force.h"

#include <algorithm>
#include <map>
#include <vector>

void NonbondedForce::init(Context& ctx) {
    build_combinded_list(ctx);
    build_charge_table(ctx);
    build_catype_table(ctx);
    init_backend(ctx);
}

void NonbondedForce::build_combinded_list(Context& ctx) {
    int n_atoms = ctx.n_atoms;
    data_.n_total = n_atoms;

    std::vector<int> atom_idx(n_atoms);
    std::vector<uint8_t> category(n_atoms);
    std::vector<int> q_state(n_atoms);
    std::vector<real_t> atom_lambdas(n_atoms);

    int cur_idx = 0;
    for (int i = 0; i < ctx.n_patoms; i++) {
        atom_idx[cur_idx] = ctx.p_atoms[i];
        category[cur_idx] = static_cast<uint8_t>(AtomCategory::P);
        q_state[cur_idx] = -1;
        atom_lambdas[cur_idx] = 1.0;
        cur_idx++;
    }

    for (int i = 0; i < ctx.n_qatoms; i++) {
        atom_idx[cur_idx] = ctx.q_atoms[i];
        category[cur_idx] = static_cast<uint8_t>(AtomCategory::Q);

        // A q-atom is dummy in a state iff its atype code is 0 (DUM is not in the
        // atypemap, so it parses to 0). In dual-topology FEP each q-atom is real in
        // exactly one state; pick that state.
        int real_state = (ctx.q_atypes[i].code != 0) ? 0 : 1;
        q_state[cur_idx] = real_state;
        atom_lambdas[cur_idx] = ctx.lambdas->cpu_data_p[real_state];

        cur_idx++;
    }

    for (int i = ctx.n_atoms_solute; i < n_atoms; i++) {
        atom_idx[cur_idx] = i;
        category[cur_idx] = static_cast<uint8_t>(AtomCategory::W);
        q_state[cur_idx] = -1;
        atom_lambdas[cur_idx] = 1.0;
        cur_idx++;
    }

    data_.atom_idx = HostDeviceBuffer<int>::from_vector(atom_idx, ctx.command_info.requested_gpu);
    data_.category = HostDeviceBuffer<uint8_t>::from_vector(category, ctx.command_info.requested_gpu);
    data_.q_state = HostDeviceBuffer<int>::from_vector(q_state, ctx.command_info.requested_gpu);
    data_.atom_lambdas = HostDeviceBuffer<real_t>::from_vector(atom_lambdas, ctx.command_info.requested_gpu);
}

void NonbondedForce::build_charge_table(Context& ctx) {
    std::vector<double> charge_table;
    charge_table.push_back(0);
    for (int i = 0; i < ctx.n_ccharges; i++) {
        charge_table.push_back(ctx.ccharges->cpu_data_p[i].charge);
    }
    for (int i = 0; i < ctx.n_qatoms; i++) {
        int real_state = (ctx.q_atypes[i].code != 0) ? 0 : 1;
        charge_table.push_back(ctx.q_charges[i + ctx.n_qatoms * real_state].charge);
    }
    std::sort(charge_table.begin(), charge_table.end());
    charge_table.resize(std::unique(charge_table.begin(), charge_table.end()) - charge_table.begin());
    std::map<double, int> charge_to_idx;
    int charge_table_size = charge_table.size();
    for (int i = 0; i < charge_table_size; i++) {
        charge_to_idx[charge_table[i]] = i;
    }

    std::vector<int> charge_types(ctx.n_atoms);

    for (int i = 0; i < ctx.n_atoms; i++) {
        int atom_idx = data_.atom_idx->cpu_data_p[i];
        auto atom_type = data_.category->cpu_data_p[i];

        if (atom_type == static_cast<uint8_t>(AtomCategory::P) ||
            atom_type == static_cast<uint8_t>(AtomCategory::W)) {
            double charge = ctx.ccharges->cpu_data_p[ctx.charges->cpu_data_p[atom_idx].code - 1].charge;
            charge_types[i] = charge_to_idx[charge];
        } else {
            int state = data_.q_state->cpu_data_p[i];
            double charge = xxx;
            charge_types[i] = charge_to_idx[charge];
        }
    }

    int sz = charge_table.size();
    std::vector<real_t> charge_pair_products(sz * sz);

    for (int i = 0; i < sz; i++) {
        for (int j = 0; j < sz; j++) {
            charge_pair_products[i * sz + j] = charge_table[i] * charge_table[j];
        }
    }

    data_.charge_types = HostDeviceBuffer<int>::from_vector(charge_types, ctx.command_info.requested_gpu);
    data_.charge_pair_products = HostDeviceBuffer<real_t>::from_vector(charge_pair_products, ctx.command_info.requested_gpu);
}

void NonbondedForce::build_catype_table(Context& ctx) {
}