#include "nonbonded_force.h"

#include <algorithm>
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
}

void NonbondedForce::build_catype_table(Context& ctx) {
}