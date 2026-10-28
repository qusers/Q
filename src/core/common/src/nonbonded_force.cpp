#include "nonbonded_force.h"

#include <algorithm>
#include <map>
#include <queue>
#include <set>
#include <utility>
#include <vector>

#include "constants.h"
#include "geometry.h"
#include "vdw_rules.h"

namespace {

struct SpatialGroupEntry {
    uint32_t morton;
    int group;
};

coord_t group_spatial_position(const Context& ctx, int group_index, uint8_t group_category) {
    const auto& config = ctx.charge_group_config;

    const auto& group = config.charge_groups[group_index];

    const coord_t* coords = ctx.coords->cpu_data_p;

    constexpr uint8_t W = static_cast<uint8_t>(AtomCategory::W);

    /*
     * Switch mode always sorts using the switch atom.
     *
     * Water also uses its switch atom when
     * iuse_switch_atom == 0, matching the CPU
     * W-W and P-W distance definitions.
     */
    if (config.iuse_switch_atom == 1 || group_category == W) {
        const int switch_atom = group.iswitch - 1;

        return coords[switch_atom];
    }

    /*
     * In all-atom mode, use the centroid for P and Q
     * groups. This is only a spatial ordering key; it
     * does not change group-pair classification.
     */
    coord_t center{0.0, 0.0, 0.0};

    for (int atom_1based : group.atoms) {
        const int atom = atom_1based - 1;

        center.x += coords[atom].x;
        center.y += coords[atom].y;
        center.z += coords[atom].z;
    }

    const double inverse_size = 1.0 / static_cast<double>(group.atoms.size());

    center.x *= inverse_size;
    center.y *= inverse_size;
    center.z *= inverse_size;

    return center;
}

void spatially_sort_group_indices(const Context& ctx, uint8_t group_category, std::vector<int>& group_indices) {
    if (group_indices.size() < 2) {
        return;
    }

    std::vector<coord_t> positions;
    positions.reserve(group_indices.size());

    for (int group : group_indices) {
        positions.push_back(group_spatial_position(ctx, group, group_category));
    }

    double minimum_x = positions[0].x;
    double minimum_y = positions[0].y;
    double minimum_z = positions[0].z;

    double maximum_x = positions[0].x;
    double maximum_y = positions[0].y;
    double maximum_z = positions[0].z;

    for (const coord_t& position : positions) {
        minimum_x = std::min(minimum_x, position.x);
        minimum_y = std::min(minimum_y, position.y);
        minimum_z = std::min(minimum_z, position.z);

        maximum_x = std::max(maximum_x, position.x);
        maximum_y = std::max(maximum_y, position.y);
        maximum_z = std::max(maximum_z, position.z);
    }

    std::vector<SpatialGroupEntry> entries;
    entries.reserve(group_indices.size());

    for (size_t i = 0; i < group_indices.size(); ++i) {
        const coord_t& position = positions[i];
        entries.push_back({get_morton_code(position, minimum_x, maximum_x, minimum_y, maximum_y, minimum_z, maximum_z), group_indices[i]});
    }

    /*
     * stable_sort keeps the original topology ordering
     * for groups with identical Morton keys.
     */
    std::stable_sort(entries.begin(), entries.end(), [](const SpatialGroupEntry& lhs, const SpatialGroupEntry& rhs) {
        return lhs.morton < rhs.morton;
    });

    for (size_t i = 0; i < entries.size(); ++i) {
        group_indices[i] = entries[i].group;
    }
}

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
    build_atom_to_group(ctx);
    build_pair_masks(ctx, parsed);
    init_backend(ctx);
}

void NonbondedForce::build_atom_to_group(Context& ctx) {
    const auto& groups = ctx.charge_group_config.charge_groups;
    const int n_groups = groups.size();

    std::vector<int> atom_to_group(ctx.n_atoms, -1);
    for (int group = 0; group < n_groups; group++) {
        for (int atom_1based : groups[group].atoms) {
            const int atom = atom_1based - 1;

            if (atom >= 0 && atom < ctx.n_atoms) {
                atom_to_group[atom] = group;
            }
        }
    }
    data_.atom_to_group = HostDeviceBuffer<int>::from_vector(atom_to_group, ctx.command_info.requested_gpu);
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
    std::vector<int> group_indices;
    std::vector<int> group_start_idx;
    std::vector<int> group_sizes;

    auto push_dummy = [&](int count) {
        for (int i = 0; i < count; i++) {
            atom_idx.push_back(-1);
            category.push_back(static_cast<uint8_t>(AtomCategory::INVALID));
            q_state.push_back(-1);
            atom_lambdas.push_back(0);
        }
    };

    std::vector<uint8_t> atom_type(ctx.n_atoms);
    for (int i = 0; i < ctx.n_patoms(); i++) {
        int idx = ctx.p_atoms[i];
        atom_type[idx] = static_cast<uint8_t>(AtomCategory::P);
    }
    for (int i = 0; i < ctx.n_qatoms(); i++) {
        int idx = ctx.q_atoms[i];
        atom_type[idx] = static_cast<uint8_t>(AtomCategory::Q);
    }
    for (int i = ctx.n_atoms_solute; i < ctx.n_atoms; i++) {
        atom_type[i] = static_cast<uint8_t>(AtomCategory::W);
    }

    const auto& groups = ctx.charge_group_config.charge_groups;
    int group_size = groups.size();

    std::vector<std::vector<int>> category_groups(3);
    /*
    Now we use switch atoms to decide the excluded atoms. Based on the three assumptions:
    1. Each atoms should be in one group.
    2. If the switch atom is exculded, all atoms that are in the same group should be excluded.
    3. The category of the atoms in one group should be same.
    */

    for (int i = 0; i < group_size; i++) {
        int atom = groups[i].iswitch - 1;
        if (ctx.excluded->cpu_data_p[atom]) continue;
        category_groups[atom_type[atom]].push_back(i);
    }

    if (ctx.command_info.requested_gpu) {
        constexpr uint8_t P = static_cast<uint8_t>(AtomCategory::P);
        constexpr uint8_t Q = static_cast<uint8_t>(AtomCategory::Q);
        constexpr uint8_t W = static_cast<uint8_t>(AtomCategory::W);

        spatially_sort_group_indices(ctx, P, category_groups[P]);
        spatially_sort_group_indices(ctx, Q, category_groups[Q]);
        spatially_sort_group_indices(ctx, W, category_groups[W]);
    }

    // P
    for (int i = 0; i < category_groups[0].size(); i++) {
        int group_idx = category_groups[0][i];
        group_indices.push_back(group_idx);
        group_start_idx.push_back(atom_idx.size());
        group_sizes.push_back(groups[group_idx].atoms.size());

        int switch_atom = groups[group_idx].iswitch - 1;
        atom_idx.push_back(switch_atom);
        category.push_back(static_cast<uint8_t>(AtomCategory::P));
        q_state.push_back(-1);
        atom_lambdas.push_back(1.0);

        for (int j = 0; j < groups[group_idx].atoms.size(); j++) {
            int atom = groups[group_idx].atoms[j] - 1;
            if (atom == switch_atom) continue;
            atom_idx.push_back(atom);
            category.push_back(static_cast<uint8_t>(AtomCategory::P));
            q_state.push_back(-1);
            atom_lambdas.push_back(1.0);
        }
    }

    int sz = atom_idx.size();
    push_dummy((32 - (sz % 32)) % 32);

    // Q
    for (int state = 0; state < ctx.n_lambdas(); state++) {
        for (int i = 0; i < category_groups[1].size(); i++) {
            int group_idx = category_groups[1][i];
            group_indices.push_back(group_idx);
            group_start_idx.push_back(atom_idx.size());
            group_sizes.push_back(groups[group_idx].atoms.size());

            int switch_atom = groups[group_idx].iswitch - 1;
            atom_idx.push_back(switch_atom);
            category.push_back(static_cast<uint8_t>(AtomCategory::Q));
            q_state.push_back(state);
            atom_lambdas.push_back(ctx.lambdas->cpu_data_p[state]);

            for (int j = 0; j < groups[group_idx].atoms.size(); j++) {
                int atom = groups[group_idx].atoms[j] - 1;
                if (atom == switch_atom) continue;
                atom_idx.push_back(atom);
                category.push_back(static_cast<uint8_t>(AtomCategory::Q));
                q_state.push_back(state);
                atom_lambdas.push_back(ctx.lambdas->cpu_data_p[state]);
            }
        }
        sz = atom_idx.size();
        push_dummy((32 - (sz % 32)) % 32);
    }

    // W
    for (int i = 0; i < category_groups[2].size(); i++) {
        int group_idx = category_groups[2][i];
        group_indices.push_back(group_idx);
        group_start_idx.push_back(atom_idx.size());
        group_sizes.push_back(groups[group_idx].atoms.size());

        int switch_atom = groups[group_idx].iswitch - 1;
        atom_idx.push_back(switch_atom);
        category.push_back(static_cast<uint8_t>(AtomCategory::W));
        q_state.push_back(-1);
        atom_lambdas.push_back(1.0);

        for (int j = 0; j < groups[group_idx].atoms.size(); j++) {
            int atom = groups[group_idx].atoms[j] - 1;
            if (atom == switch_atom) continue;
            atom_idx.push_back(atom);
            category.push_back(static_cast<uint8_t>(AtomCategory::W));
            q_state.push_back(-1);
            atom_lambdas.push_back(1.0);
        }
    }
    sz = atom_idx.size();
    push_dummy((32 - (sz % 32)) % 32);

    sz = atom_idx.size();
    data_.n_total = sz;

    data_.atom_idx = HostDeviceBuffer<int>::from_vector(atom_idx, ctx.command_info.requested_gpu);
    data_.category = HostDeviceBuffer<uint8_t>::from_vector(category, ctx.command_info.requested_gpu);
    data_.q_state = HostDeviceBuffer<int>::from_vector(q_state, ctx.command_info.requested_gpu);
    data_.atom_lambdas = HostDeviceBuffer<real_t>::from_vector(atom_lambdas, ctx.command_info.requested_gpu);
    data_.group_indices = HostDeviceBuffer<int>::from_vector(group_indices, ctx.command_info.requested_gpu);
    data_.group_start_idx = HostDeviceBuffer<int>::from_vector(group_start_idx, ctx.command_info.requested_gpu);
    data_.group_sizes = HostDeviceBuffer<int>::from_vector(group_sizes, ctx.command_info.requested_gpu);
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
