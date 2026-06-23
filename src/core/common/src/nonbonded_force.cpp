#include "nonbonded_force.h"

#include <algorithm>
#include <map>
#include <vector>

#include "constants.h"
#include "vdw_rules.h"

void NonbondedForce::init(Context& ctx) {
    build_combinded_list(ctx);
    build_charge_table(ctx);
    build_catype_table(ctx);
    init_backend(ctx);
}

void NonbondedForce::build_combinded_list(Context& ctx) {
    std::vector<int> atom_idx;
    std::vector<uint8_t> category;
    std::vector<int> q_state;
    std::vector<real_t> atom_lambdas;

    auto push_dummy = [&](int count) {
        for (int i = 0; i < count; i++) {
            atom_idx.push_back(-1);
            category.push_back(static_cast<uint8_t>(AtomCategory::INVALID));
            q_state.push_back(-1);
            atom_lambdas.push_back(0);
        }
    };

    for (int i = 0; i < ctx.n_patoms; i++) {
        int idx = ctx.p_atoms[i];
        if (ctx.excluded->cpu_data_p[idx]) continue;
        atom_idx.push_back(idx);
        category.push_back(static_cast<uint8_t>(AtomCategory::P));
        q_state.push_back(-1);
        atom_lambdas.push_back(1.0);
    }

    int sz = atom_idx.size();
    push_dummy((32 - (sz % 32)) % 32);

    for (int state = 0; state < ctx.n_lambdas; state++) {
        for (int i = 0; i < ctx.n_qatoms; i++) {
            int idx = ctx.q_atoms[i];
            if (ctx.excluded->cpu_data_p[idx]) continue;
            atom_idx.push_back(idx);
            category.push_back(static_cast<uint8_t>(AtomCategory::Q));
            q_state.push_back(state);
            atom_lambdas.push_back(ctx.lambdas->cpu_data_p[state]);
        }
        sz = atom_idx.size();
        push_dummy((32 - (sz % 32)) % 32);
    }

    for (int i = ctx.n_atoms_solute; i < ctx.n_atoms; i++) {
        if (ctx.excluded->cpu_data_p[i]) continue;
        atom_idx.push_back(i);
        category.push_back(static_cast<uint8_t>(AtomCategory::W));
        q_state.push_back(-1);
        atom_lambdas.push_back(1.0);
    }

    sz = atom_idx.size();
    push_dummy((32 - (sz % 32)) % 32);

    sz = atom_idx.size();
    data_.n_total = sz;
    data_.atom_idx = HostDeviceBuffer<int>::from_vector(atom_idx, ctx.command_info.requested_gpu);
    data_.category = HostDeviceBuffer<uint8_t>::from_vector(category, ctx.command_info.requested_gpu);
    data_.q_state = HostDeviceBuffer<int>::from_vector(q_state, ctx.command_info.requested_gpu);
    data_.atom_lambdas = HostDeviceBuffer<real_t>::from_vector(atom_lambdas, ctx.command_info.requested_gpu);
}

void NonbondedForce::build_charge_table(Context& ctx) {
    std::map<int, int> atom_idx_to_q_idx;
    for (int i = 0; i < ctx.n_qatoms; i++) {
        int atom_idx = ctx.q_atoms[i];
        atom_idx_to_q_idx[atom_idx] = i;
    }

    std::vector<real_t> atom_charge(data_.n_total);
    for (int i = 0; i < data_.n_total; i++) {
        int atom_idx = data_.atom_idx->cpu_data_p[i];
        auto atom_type = data_.category->cpu_data_p[i];

        if (atom_type == static_cast<uint8_t>(AtomCategory::INVALID)) {
            atom_charge[i] = 0;
            continue;
        }

        if (atom_type == static_cast<uint8_t>(AtomCategory::P) ||
            atom_type == static_cast<uint8_t>(AtomCategory::W)) {
            double charge = ctx.ccharges->cpu_data_p[ctx.charges->cpu_data_p[atom_idx].code - 1].charge;
            atom_charge[i] = charge;
        } else {
            int state = data_.q_state->cpu_data_p[i];
            int q_idx = atom_idx_to_q_idx[atom_idx];
            double charge = ctx.q_charges[q_idx + ctx.n_qatoms * state].charge;
            atom_charge[i] = charge;
        }
    }
    data_.atom_charge = HostDeviceBuffer<real_t>::from_vector(atom_charge, ctx.command_info.requested_gpu);
}

void NonbondedForce::build_catype_table(Context& ctx) {
    std::vector<vdw_atom_param_t> atom_vdw(data_.n_total);
    auto& catypes = ctx.catypes->cpu_data_p;

    std::map<int, int> atom_idx_to_q_idx;
    for (int i = 0; i < ctx.n_qatoms; i++) {
        int atom_idx = ctx.q_atoms[i];
        atom_idx_to_q_idx[atom_idx] = i;
    }

    for (int i = 0; i < data_.n_total; i++) {
        int atom_idx = data_.atom_idx->cpu_data_p[i];
        auto atom_type = data_.category->cpu_data_p[i];

        if (atom_type == static_cast<uint8_t>(AtomCategory::INVALID)) {
            atom_vdw[i] = vdw_atom_param_t{0, 0, 0, 0};
            continue;
        }

        if (atom_type == static_cast<uint8_t>(AtomCategory::P) || atom_type == static_cast<uint8_t>(AtomCategory::W)) {
            const catype_t& catype = catypes[ctx.atypes->cpu_data_p[atom_idx].code - 1];
            atom_vdw[i] = vdw_atom_param_t{catype.aii_normal, catype.bii_normal, catype.aii_1_4, catype.bii_1_4};
        } else {
            int state = data_.q_state->cpu_data_p[i];
            int q_idx = atom_idx_to_q_idx[atom_idx];
            const atype_t& atype = ctx.q_atypes[q_idx + ctx.n_qatoms * state];
            if (atype.code > 0) {
                const catype_t& catype = ctx.q_catypes[atype.code - 1];
                atom_vdw[i] = vdw_atom_param_t{catype.aii_normal, catype.bii_normal, catype.aii_1_4, catype.bii_1_4};
            } else {
                catype_t zero = {};
                atom_vdw[i] = vdw_atom_param_t{0, 0, 0, 0};
            }
        }
    }
    data_.atom_vdw = HostDeviceBuffer<vdw_atom_param_t>::from_vector(atom_vdw, ctx.command_info.requested_gpu);
}
