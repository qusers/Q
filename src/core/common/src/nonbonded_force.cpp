#include "nonbonded_force.h"

#include <algorithm>
#include <map>
#include <queue>
#include <set>
#include <vector>

#include "constants.h"
#include "vdw_rules.h"

namespace {

typedef std::array<int, 3> Status;

std::pair<std::set<Status>, std::set<Status>> get_solute_excluded_and_34_pairs(Context& ctx, const ParseResult& parsed) {
    const int nsolute = ctx.n_atoms_solute;
    const int n_lambda = ctx.n_lambdas();

    std::vector<int> atom_to_q(nsolute, -1);
    const int nqatoms = ctx.n_qatoms();
    for (int q = 0; q < nqatoms; q++) {
        atom_to_q[ctx.q_atoms[q]] = q;
    }

    const size_t q_bonds_per_state = parsed.q_bonds.size() / n_lambda;
    std::set<Status> result_excluded, result_34;

    for (int s = 0; s < n_lambda; s++) {
        std::vector<std::vector<int>> g(nsolute);
        std::set<std::pair<int, int>> q_bonds;

        int start = q_bonds_per_state * s;
        for (int i = start; i < start + q_bonds_per_state; i++) {
            const auto& bond = parsed.q_bonds[i];
            int ai = bond.ai - 1;
            int aj = bond.aj - 1;
            if (ai > aj) {
                std::swap(ai, aj);
            }

            if (bond.code > 0) {
                g[ai].push_back(aj);
                g[aj].push_back(ai);
            }
            q_bonds.insert({ai, aj});
        }

        for (int i = 0; i < parsed.n_bonds_solute; i++) {
            const auto& bond = parsed.bonds[i];
            if (bond.code == 0) continue;
            int ai = bond.ai - 1, aj = bond.aj - 1;
            if (ai > aj) {
                std::swap(ai, aj);
            }
            if (q_bonds.count({ai, aj})) continue;

            g[ai].push_back(aj);
            g[aj].push_back(ai);
        }

        std::queue<int> Q;
        std::vector<int> distance(nsolute, -1);

        for (const int origin : ctx.q_atoms) {
            std::fill(distance.begin(), distance.end(), -1);
            distance[origin] = 0;
            Q.push(origin);

            while (!Q.empty()) {
                int u = Q.front();
                Q.pop();
                if (distance[u] == 3) continue;

                for (int v : g[u]) {
                    if (distance[v] >= 0) continue;

                    const int dis = distance[u] + 1;
                    distance[v] = dis;
                    Q.push(v);

                    const int ai = std::min(origin, v);
                    const int aj = std::max(origin, v);
                    if (dis <= 2) {
                        result_excluded.insert({ai, aj, s});
                    } else {
                        result_34.insert({ai, aj, s});
                    }
                }
            }
        }
    }

    // add q_exclpairs

    const size_t exclusion_per_state = parsed.q_exclpairs.size() / n_lambda;
    for (int i = 0; i < parsed.q_exclpairs.size(); i++) {
        int state = i / exclusion_per_state;
        int ai = parsed.q_exclpairs[i].ai - 1, aj = parsed.q_exclpairs[i].aj - 1, excl = parsed.q_exclpairs[i].excl;
        if (ai > aj) std::swap(ai, aj);

        if (excl) {
            result_excluded.insert({ai, aj, state});
        }
    }

    // For p-p, it should be same as cur_g
    // Do one check in here

    // Physical solute atom index -> Q-atom index.

    auto p_type = static_cast<uint8_t>(AtomCategory::P);
    for (int i = 0; i < nsolute; i++) {
        if (atom_to_q[i] >= 0) continue;
        for (int j = i + 1; j < nsolute; j++) {
            if (atom_to_q[j] >= 0) continue;
            auto bond_type = get_bond_type(nsolute, ctx.LJ_matrix->cpu_data_p, i, p_type, j, p_type);

            if (bond_type == BondType::Bond23) {
                result_excluded.insert({i, j, 0});
            } else if (bond_type == BondType::Bond14) {
                result_34.insert({i, j, 0});
            }
        }
    }

    return {result_excluded, result_34};
}

}  // namespace

void NonbondedForce::init(Context& ctx, const ParseResult& parsed) {
    build_combinded_list(ctx);
    build_charge_table(ctx);
    build_catype_table(ctx);
    build_pair_masks(ctx, parsed);
    init_backend(ctx);
}

void NonbondedForce::build_pair_masks(Context& ctx, const ParseResult& parsed) {
    auto* atom_idx = data_.atom_idx->cpu_data_p;
    auto* category = data_.category->cpu_data_p;
    auto* q_state = data_.q_state->cpu_data_p;

    const bool run_gpu = ctx.command_info.requested_gpu;
    int nsolute = parsed.n_atoms_solute;

    auto [result_excluded, result_34] = get_solute_excluded_and_34_pairs(ctx, parsed);
    data_.n_cluster = data_.n_total / 32;
    const int n = data_.n_cluster;
    data_.n_pair_tiles = n * (n + 1) / 2;

    std::vector<uint32_t> pair_masks(data_.n_pair_tiles * 32, 0);
    std::vector<uint32_t> pair_14_masks(data_.n_pair_tiles * 32, 0);
    auto Q = static_cast<uint8_t>(AtomCategory::Q);
    auto W = static_cast<uint8_t>(AtomCategory::W);
    auto P = static_cast<uint8_t>(AtomCategory::P);

    auto classify_pair = [&](int x_slot, int y_slot) -> CalculationType {
        int ai = atom_idx[x_slot];
        int aj = atom_idx[y_slot];

        if (ai < 0 || aj < 0 || ai == aj) {
            return CalculationType::Exclude;
        }

        if (ai > aj) {
            std::swap(ai, aj);
            std::swap(x_slot, y_slot);
        }

        const uint8_t ci = category[x_slot];
        const bool qi = ci == Q;
        const bool wi = ci == W;
        const bool pi = ci == P;

        const uint8_t cj = category[y_slot];
        const bool qj = cj == Q;
        const bool wj = cj == W;
        const bool pj = cj == P;

        if (qi && qj && q_state[x_slot] != q_state[y_slot]) {
            return CalculationType::Exclude;
        }

        if (wi || wj) {
            if (wi && wj && (ai - nsolute) / 3 == (aj - nsolute) / 3) {
                return CalculationType::Exclude;
            }
            return CalculationType::Calculate;
        }

        int state;
        if (qi && qj) {
            // q q
            state = q_state[x_slot];
        } else if (pi && pj) {
            // p p
            state = 0;
        } else {
            // p q
            if (qi) {
                state = q_state[x_slot];
            } else {
                state = q_state[y_slot];
            }
        }

        if (result_excluded.count({ai, aj, state})) {
            return CalculationType::Exclude;
        } else if (result_34.count({ai, aj, state})) {
            return CalculationType::Calculate_14;
        } else {
            return CalculationType::Calculate;
        }
    };

    for (int tx = 0; tx < data_.n_cluster; tx++) {
        for (int ty = tx; ty < data_.n_cluster; ty++) {
            const int tile = pair_tile_index(tx, ty, data_.n_cluster);

            for (int x_lane = 0; x_lane < 32; x_lane++) {
                const int x_slot = tx * 32 + x_lane;

                uint32_t active = 0;
                uint32_t is14 = 0;

                for (int y_lane = 0; y_lane < 32; y_lane++) {
                    if (tx == ty && y_lane <= x_lane) continue;

                    const int y_slot = ty * 32 + y_lane;
                    auto calculation_type = classify_pair(x_slot, y_slot);

                    if (calculation_type == CalculationType::Exclude) continue;

                    const uint32_t bit = uint32_t{1} << y_lane;

                    active |= bit;

                    if (calculation_type == CalculationType::Calculate_14) {
                        is14 |= bit;
                    }
                }
                const size_t row = tile * 32 + x_lane;
                pair_masks[row] = active;
                pair_14_masks[row] = is14;
            }
        }
    }

    data_.pair_masks = HostDeviceBuffer<uint32_t>::from_vector(pair_masks, run_gpu);
    data_.pair_14_masks = HostDeviceBuffer<uint32_t>::from_vector(pair_14_masks, run_gpu);
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

    for (int i = 0; i < ctx.n_patoms(); i++) {
        int idx = ctx.p_atoms[i];
        if (ctx.excluded->cpu_data_p[idx]) continue;
        atom_idx.push_back(idx);
        category.push_back(static_cast<uint8_t>(AtomCategory::P));
        q_state.push_back(-1);
        atom_lambdas.push_back(1.0);
    }

    int sz = atom_idx.size();
    push_dummy((32 - (sz % 32)) % 32);

    for (int state = 0; state < ctx.n_lambdas(); state++) {
        for (int i = 0; i < ctx.n_qatoms(); i++) {
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
    for (int i = 0; i < ctx.n_qatoms(); i++) {
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
            double charge = ctx.q_charges[q_idx + ctx.n_qatoms() * state].charge;
            atom_charge[i] = charge;
        }
    }
    data_.atom_charge = HostDeviceBuffer<real_t>::from_vector(atom_charge, ctx.command_info.requested_gpu);
}

void NonbondedForce::build_catype_table(Context& ctx) {
    std::vector<vdw_atom_param_t> atom_vdw(data_.n_total);
    auto& catypes = ctx.catypes->cpu_data_p;

    std::map<int, int> atom_idx_to_q_idx;
    for (int i = 0; i < ctx.n_qatoms(); i++) {
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
            const atype_t& atype = ctx.q_atypes[q_idx + ctx.n_qatoms() * state];
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
