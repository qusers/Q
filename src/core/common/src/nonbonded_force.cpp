#include "nonbonded_force.h"

#include <algorithm>
#include <map>
#include <queue>
#include <set>
#include <utility>
#include <vector>

#include "constants.h"
#include "vdw_rules.h"

namespace {

typedef std::array<int, 3> Status;

struct SolutePairClassification {
    // Pairs skipped by the ordinary NB kernel
    // physical exclusions + migrated special pairs;
    std::set<Status> excluded;
    std::set<Status> pairs_14;
    std::set<Status> special;
};

SolutePairClassification classify_solute_pairs(Context& ctx, const ParseResult& parsed) {
    const int nsolute = ctx.n_atoms_solute;
    const int n_lambda = ctx.n_lambdas();

    std::vector<int> atom_to_q(nsolute, -1);
    const int nqatoms = ctx.n_qatoms();
    for (int q = 0; q < nqatoms; q++) {
        atom_to_q[ctx.q_atoms[q]] = q;
    }

    const size_t q_bonds_per_state = parsed.q_bonds.size() / n_lambda;
    std::set<Status> result_excluded, result_14, result_special;

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
                        result_14.insert({ai, aj, s});
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
                result_14.insert({i, j, 0});
            }
        }
    }

    // Find special pairs, currently for softpair
    for (const auto& pair : parsed.q_softpairs) {
        int ai = ctx.q_atoms.at(pair.qi - 1);
        int aj = ctx.q_atoms.at(pair.qj - 1);

        if (ai > aj) std::swap(ai, aj);

        // Match the atoms omitted from the combined NB list.
        if (ctx.excluded->cpu_data_p[ai] || ctx.excluded->cpu_data_p[aj]) {
            continue;
        }

        for (int state = 0; state < n_lambda; ++state) {
            const Status key{ai, aj, state};

            // Preserve physical exclusions and 1–4 priority.
            if (result_excluded.count(key) || result_14.count(key)) {
                continue;
            }

            result_special.insert(key);
            result_excluded.insert(key);
        }
    }

    return {std::move(result_excluded), std::move(result_14), std::move(result_special)};
}

std::vector<special_pair_term_t> build_lj14_pairs(Context& ctx, const std::set<Status>& pairs_14, const std::set<Status>& excluded) {
    const int nqat = ctx.n_qatoms();
    std::vector<int> atom_to_q(ctx.n_atoms_solute, -1);

    for (int q = 0; q < nqat; ++q) {
        atom_to_q.at(ctx.q_atoms[q]) = q;
    }
    auto atom_params = [&](int atom, int state) -> std::pair<double, const catype_t*> {
        const int q = atom_to_q.at(atom);

        if (q >= 0) {
            const size_t index = static_cast<size_t>(state) * nqat + q;
            const int type = ctx.q_atypes.at(index).code;

            return {
                ctx.q_charges.at(index).charge,
                &ctx.q_catypes.at(type - 1)};
        }

        const int charge_code = ctx.charges->cpu_data_p[atom].code;
        const int type_code = ctx.atypes->cpu_data_p[atom].code;

        return {
            ctx.ccharges->cpu_data_p[charge_code - 1].charge,
            &ctx.catypes->cpu_data_p[type_code - 1]};
    };

    std::vector<special_pair_term_t> terms;
    terms.reserve(pairs_14.size());

    for (const auto& key : pairs_14) {
        // An explicit exclusion overrides 1–4 classification.
        if (excluded.count(key)) continue;

        const int ai = key[0];
        const int aj = key[1];
        const int s = key[2];

        if (ctx.excluded->cpu_data_p[ai] || ctx.excluded->cpu_data_p[aj]) {
            continue;
        }

        const bool qi = atom_to_q.at(ai) >= 0;
        const bool qj = atom_to_q.at(aj) >= 0;

        const auto [charge_i, param_i] = atom_params(ai, s);
        const auto [charge_j, param_j] = atom_params(aj, s);

        const auto pair = combine_vdw(ctx.topo.vdw_rule, param_i->aii_1_4, param_i->bii_1_4, param_j->aii_1_4, param_j->bii_1_4);

        special_pair_term_t term{};
        term.ai = ai;
        term.aj = aj;
        term.state = (qi || qj) ? s : -1;
        term.kind = SpecialPairKind::LJ14;

        const auto category_i = static_cast<uint8_t>(qi ? AtomCategory::Q : AtomCategory::P);
        const auto category_j = static_cast<uint8_t>(qj ? AtomCategory::Q : AtomCategory::P);

        term.coul_slot = nb_coul_slot(category_i, category_j, qi ? s : -1, qj ? s : -1, ctx.n_lambdas());
        term.vdw_slot = term.coul_slot + 1;

        term.coulomb_coeff = ctx.topo.coulomb_constant * ctx.topo.el14_scale * charge_i * charge_j;

        term.p1 = pair.x;
        term.p2 = pair.y;

        terms.push_back(term);
    }

    return terms;
}

std::vector<special_pair_term_t> build_special_pairs(Context& ctx, const ParseResult& parsed, const std::set<Status>& s) {
    std::vector<special_pair_term_t> res;
    res.reserve(s.size());
    int nsolute = parsed.n_atoms_solute;
    std::vector<int> atom_to_q(nsolute, -1);
    const int nqat = ctx.n_qatoms();

    for (int q = 0; q < nqat; q++) {
        atom_to_q.at(ctx.q_atoms[q]) = q;
    }

    for (const auto& key : s) {
        const int ai = key[0];
        const int aj = key[1];
        const int state = key[2];

        const int iq = atom_to_q.at(ai);
        const int jq = atom_to_q.at(aj);

        if (iq < 0 || jq < 0) {
            throw std::runtime_error("Soft pairs must contain two Q-atoms.");
        }

        const size_t index_i = static_cast<size_t>(state) * nqat + iq;
        const size_t index_j = static_cast<size_t>(state) * nqat + jq;

        const int type_i = ctx.q_atypes.at(index_i).code;
        const int type_j = ctx.q_atypes.at(index_j).code;

        if (type_i <= 0 || type_j <= 0) {
            throw std::runtime_error("Soft pairs require defined Q-atom types.");
        }

        const auto& param_i = ctx.q_catypes.at(type_i - 1);
        const auto& param_j = ctx.q_catypes.at(type_j - 1);

        const double charge_i = ctx.q_charges.at(index_i).charge;
        const double charge_j = ctx.q_charges.at(index_j).charge;

        special_pair_term_t term{};
        term.ai = ai;
        term.aj = aj;
        term.state = state;
        term.kind = SpecialPairKind::SoftRepulsion;

        term.coul_slot = EnergyBuffer::eq_index(ENERGY_FIXED_COUNT, state, EQ_NB_QQ_COUL);
        term.vdw_slot = EnergyBuffer::eq_index(ENERGY_FIXED_COUNT, state, EQ_NB_QQ_VDW);

        term.coulomb_coeff = ctx.topo.coulomb_constant * charge_i * charge_j;

        term.p1 = param_i.soft_C * param_j.soft_C;
        term.p2 = param_i.soft_a * param_j.soft_a;

        res.emplace_back(term);
    }
    return std::move(res);
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

    auto [result_excluded, result_14, result_special] = classify_solute_pairs(ctx, parsed);

    data_.special_pairs = build_special_pairs(ctx, parsed, result_special);
    const auto pairs14 = build_lj14_pairs(ctx, result_14, result_excluded);
    data_.special_pairs.insert(data_.special_pairs.end(), pairs14.begin(), pairs14.end());

    data_.n_cluster = data_.n_total / 32;
    const int n = data_.n_cluster;
    data_.n_pair_tiles = n * (n + 1) / 2;

    std::vector<uint32_t> pair_masks(data_.n_pair_tiles * 32, 0);
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
        } else if (result_14.count({ai, aj, state})) {
            return CalculationType::Exclude;
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

                for (int y_lane = 0; y_lane < 32; y_lane++) {
                    if (tx == ty && y_lane <= x_lane) continue;

                    const int y_slot = ty * 32 + y_lane;
                    auto calculation_type = classify_pair(x_slot, y_slot);

                    if (calculation_type == CalculationType::Exclude) continue;

                    const uint32_t bit = uint32_t{1} << y_lane;

                    active |= bit;
                }
                const size_t row = tile * 32 + x_lane;
                pair_masks[row] = active;
            }
        }
    }

    data_.pair_masks = HostDeviceBuffer<uint32_t>::from_vector(pair_masks, run_gpu);
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
            atom_vdw[i] = vdw_atom_param_t{0, 0};
            continue;
        }

        if (atom_type == static_cast<uint8_t>(AtomCategory::P) || atom_type == static_cast<uint8_t>(AtomCategory::W)) {
            const catype_t& catype = catypes[ctx.atypes->cpu_data_p[atom_idx].code - 1];
            atom_vdw[i] = vdw_atom_param_t{catype.aii_normal, catype.bii_normal};
        } else {
            int state = data_.q_state->cpu_data_p[i];
            int q_idx = atom_idx_to_q_idx[atom_idx];
            const atype_t& atype = ctx.q_atypes[q_idx + ctx.n_qatoms() * state];
            if (atype.code > 0) {
                const catype_t& catype = ctx.q_catypes[atype.code - 1];
                atom_vdw[i] = vdw_atom_param_t{catype.aii_normal, catype.bii_normal};
            } else {
                catype_t zero = {};
                atom_vdw[i] = vdw_atom_param_t{0, 0};
            }
        }
    }
    data_.atom_vdw = HostDeviceBuffer<vdw_atom_param_t>::from_vector(atom_vdw, ctx.command_info.requested_gpu);
}
