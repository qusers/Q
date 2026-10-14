#include "nonbonded_force.h"

#include <algorithm>
#include <map>
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
}  // namespace

void NonbondedForce::init(Context& ctx) {
    build_combinded_list(ctx);
    build_charge_table(ctx);
    build_catype_table(ctx);
    build_atom_to_group(ctx);
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
