#include "nonbonded_force.h"

#include <algorithm>
#include <map>
#include <vector>

#include "constants.h"
#include "vdw_rules.h"

namespace {
typedef std::array<real_t, 4> catype_key;
catype_key build_catype_key(const catype_t& catype) {
    return {catype.aii_normal, catype.bii_normal, catype.aii_1_4, catype.bii_1_4};
}

int add_catype(const catype_t& catype, std::map<catype_key, int>& mp, std::vector<catype_t>& table) {
    const auto& key = build_catype_key(catype);
    if (mp.find(key) == mp.end()) {
        int sz = static_cast<int>(table.size());
        table.push_back(catype);
        mp[key] = sz;
    }
    return mp[key];
}

}  // namespace

void NonbondedForce::init(Context& ctx) {
    build_combinded_list(ctx);
    build_charge_table(ctx);
    build_catype_table(ctx);
    init_backend(ctx);
}

void NonbondedForce::build_combinded_list(Context& ctx) {
    int n_atoms = ctx.n_atoms;
    data_.n_total = n_atoms - ctx.n_excluded;

    std::vector<int> atom_idx(data_.n_total);
    std::vector<uint8_t> category(data_.n_total);
    std::vector<int> q_state(data_.n_total);
    std::vector<real_t> atom_lambdas(data_.n_total);

    int cur_idx = 0;
    for (int i = 0; i < ctx.n_patoms; i++) {
        int idx = ctx.p_atoms[i];
        if (ctx.excluded->cpu_data_p[idx]) continue;
        atom_idx[cur_idx] = idx;
        category[cur_idx] = static_cast<uint8_t>(AtomCategory::P);
        q_state[cur_idx] = -1;
        atom_lambdas[cur_idx] = 1.0;
        cur_idx++;
    }

    for (int i = 0; i < ctx.n_qatoms; i++) {
        int idx = ctx.q_atoms[i];
        if (ctx.excluded->cpu_data_p[idx]) continue;
        atom_idx[cur_idx] = idx;
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
        if (ctx.excluded->cpu_data_p[i]) continue;
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

    std::vector<int> charge_types(data_.n_total);

    std::map<int, int> atom_idx_to_q_idx;
    for (int i = 0; i < ctx.n_qatoms; i++) {
        int atom_idx = ctx.q_atoms[i];
        atom_idx_to_q_idx[atom_idx] = i;
    }

    for (int i = 0; i < data_.n_total; i++) {
        int atom_idx = data_.atom_idx->cpu_data_p[i];
        auto atom_type = data_.category->cpu_data_p[i];

        if (atom_type == static_cast<uint8_t>(AtomCategory::P) ||
            atom_type == static_cast<uint8_t>(AtomCategory::W)) {
            double charge = ctx.ccharges->cpu_data_p[ctx.charges->cpu_data_p[atom_idx].code - 1].charge;
            charge_types[i] = charge_to_idx[charge];
        } else {
            int state = data_.q_state->cpu_data_p[i];
            int q_idx = atom_idx_to_q_idx[atom_idx];
            double charge = ctx.q_charges[q_idx + ctx.n_qatoms * state].charge;
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
    data_.zero_charge_type = charge_to_idx[0];
}

void NonbondedForce::build_catype_table(Context& ctx) {
    std::vector<catype_t> all_catypes;
    std::map<catype_key, int> catype_key_to_idx;

    auto& catypes = ctx.catypes->cpu_data_p;

    for (int i = 0; i < ctx.n_catypes; i++) {
        const auto& catype = catypes[i];
        add_catype(catype, catype_key_to_idx, all_catypes);
    }

    for (int i = 0; i < ctx.n_qatoms; i++) {
        int real_state = (ctx.q_atypes[i].code != 0) ? 0 : 1;
        const auto& atype = ctx.q_atypes[i + ctx.n_qatoms * real_state];
        const auto& catype = ctx.q_catypes[atype.code - 1];
        add_catype(catype, catype_key_to_idx, all_catypes);
    }
    // Put zero_catype in it
    catype_t zero_catype = {};
    add_catype(zero_catype, catype_key_to_idx, all_catypes);

    std::map<int, int> atom_idx_to_q_idx;
    for (int i = 0; i < ctx.n_qatoms; i++) {
        int atom_idx = ctx.q_atoms[i];
        atom_idx_to_q_idx[atom_idx] = i;
    }

    std::vector<int> catype_types(data_.n_total);

    for (int i = 0; i < data_.n_total; i++) {
        int atom_idx = data_.atom_idx->cpu_data_p[i];
        auto atom_type = data_.category->cpu_data_p[i];
        if (atom_type == static_cast<uint8_t>(AtomCategory::P) || atom_type == static_cast<uint8_t>(AtomCategory::W)) {
            const catype_t& catype = catypes[ctx.atypes->cpu_data_p[atom_idx].code - 1];
            catype_types[i] = catype_key_to_idx[build_catype_key(catype)];
        } else {
            int state = data_.q_state->cpu_data_p[i];
            int q_idx = atom_idx_to_q_idx[atom_idx];
            const atype_t& atype = ctx.q_atypes[q_idx + ctx.n_qatoms * state];
            const catype_t& catype = ctx.q_catypes[atype.code - 1];
            catype_types[i] = catype_key_to_idx[build_catype_key(catype)];
        }
    }

    int sz = all_catypes.size();
    std::vector<vdw_pair_param_t> catype_pair_params(sz * sz);

    for (int i = 0; i < sz; i++) {
        for (int j = 0; j < sz; j++) {
            const catype_t& ci = all_catypes[i];
            const catype_t& cj = all_catypes[j];
            vdw_pair_param_t pair_param = {};
            if (ctx.topo.vdw_rule == VDW_GEOMETRIC) {
                calc_vdw_geometric(ci.aii_normal, cj.aii_normal, ci.bii_normal, cj.bii_normal, static_cast<real_t>(1.0), &pair_param.a_normal, &pair_param.b_normal);
                calc_vdw_geometric(ci.aii_1_4, cj.aii_1_4, ci.bii_1_4, cj.bii_1_4, static_cast<real_t>(1.0), &pair_param.a_14, &pair_param.b_14);
            } else {
                calc_vdw_arithmetic(ci.aii_normal, cj.aii_normal, ci.bii_normal, cj.bii_normal, static_cast<real_t>(1.0), &pair_param.a_normal, &pair_param.b_normal);
                calc_vdw_arithmetic(ci.aii_1_4, cj.aii_1_4, ci.bii_1_4, cj.bii_1_4, static_cast<real_t>(1.0), &pair_param.a_14, &pair_param.b_14);
            }
            catype_pair_params[i * sz + j] = pair_param;
        }
    }

    data_.catype_types = HostDeviceBuffer<int>::from_vector(catype_types, ctx.command_info.requested_gpu);
    data_.catype_pair_params = HostDeviceBuffer<vdw_pair_param_t>::from_vector(catype_pair_params, ctx.command_info.requested_gpu);
    data_.zero_catype_type = catype_key_to_idx[build_catype_key(zero_catype)];
}
