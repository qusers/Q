#include "cpu_nonbonded_force.h"

#include <vector>

#include "constants.h"
#include "cpu_force_accumulation.h"
#include "geometry.h"

namespace {
void accumulate_energy(Context& ctx, real_t vel, real_t vvdw,
                       uint8_t atom1_type, uint8_t atom2_type, int atom1_state, int atom2_state) {
    energy_accum_t* e = ctx.energy.host();
    int coul = nb_coul_slot(atom1_type, atom2_type, atom1_state, atom2_state, ctx.n_lambdas());
    add_energy(e[coul], vel);
    add_energy(e[coul + 1], vvdw);  // vdw slot is adjacent (same invariant as GPU)
}

void accumulate_lrf_source(LrfCoefficients& target, const coord_t& source_coord, double source_charge) {
    coord_t r = source_coord - target.center;
    double r_len2 = norm2(r);
    double r_len = std::sqrt(r_len2);
    double inv_r_len = 1.0 / r_len;
    double inv_r_len2 = 1.0 / r_len2;
    double inv_r_len3 = inv_r_len * inv_r_len2;
    double inv_r_len5 = inv_r_len3 * inv_r_len2;
    double inv_r_len7 = inv_r_len5 * inv_r_len2;

    target.phi0 += source_charge * inv_r_len;  // q / r

    double r_array[3] = {r.x, r.y, r.z};
    for (int i = 0; i < 3; i++) {
        target.phi1[i] -= source_charge * r_array[i] * inv_r_len3;
    }

    for (int i = 0; i < 3; i++) {
        for (int j = 0; j < 3; j++) {
            double delta = i == j ? inv_r_len3 : 0;
            target.phi2[i * 3 + j] += source_charge * ((3.0 * r_array[i] * r_array[j] * inv_r_len5) - delta);
        }
    }

    for (int a = 0; a < 3; a++) {
        for (int b = 0; b < 3; b++) {
            for (int c = 0; c < 3; c++) {
                int idx = (a * 3 + b) * 3 + c;

                double delta_ab = a == b;
                double delta_ac = a == c;
                double delta_bc = b == c;

                double v1 = 3.0 * (delta_ab * r_array[c] + delta_ac * r_array[b] + delta_bc * r_array[a]) * inv_r_len5;
                double v2 = -15.0 * r_array[a] * r_array[b] * r_array[c] * inv_r_len7;

                target.phi3[idx] += source_charge * (v1 + v2);
            }
        }
    }
}

}  // namespace

void CpuNonbondedForce::calc_all_direct_pairs(Context& ctx) {
    int sz = data_.n_total;
    for (int i = 0; i < sz; i++) {
        for (int j = i + 1; j < sz; j++) {
            calc_direct_pair(ctx, i, j);
        }
    }
}

void CpuNonbondedForce::init_calculation_groups(Context& ctx) {
    const auto& config = ctx.charge_group_config;
    const auto& groups = config.charge_groups;

    const int n_groups = groups.size();
    const int n_solute_groups = config.n_cgrps_solute;

    const coord_t* coords = ctx.coords->cpu_data_p;
    const bool* excluded = ctx.excluded->cpu_data_p;

    auto normal_cutoff = [&](int group1, int group2) {
        const bool solute1 = group1 < n_solute_groups;
        const bool solute2 = group2 < n_solute_groups;
        if (solute1 && solute2) {
            return ctx.md.solute_solute;
        } else if (!solute1 && !solute2) {
            return ctx.md.solvent_solvent;
        } else {
            return ctx.md.solute_solvent;
        }
    };

    auto group_distance2 = [&](int group1, int group2) {
        if (group1 == group2) return 0.0;
        const bool solute1 = group1 < n_solute_groups;
        const bool solute2 = group2 < n_solute_groups;

        if (config.iuse_switch_atom == 1) {
            const int atom1 = groups[group1].iswitch - 1;
            const int atom2 = groups[group2].iswitch - 1;
            return norm2(coords[atom1] - coords[atom2]);
        }

        if (!solute1 && !solute2) {
            const int atom1 = groups[group1].iswitch - 1;
            const int atom2 = groups[group2].iswitch - 1;
            return norm2(coords[atom1] - coords[atom2]);
        }

        if (solute1 != solute2) {
            const int solute_group = solute1 ? group1 : group2;
            const int water_group = solute1 ? group2 : group1;

            const int water_switch = groups[water_group].iswitch - 1;

            double mi = std::numeric_limits<double>::infinity();

            for (int atom : groups[solute_group].atoms) {
                const int atom_idx = atom - 1;
                mi = std::min(mi, norm2(coords[atom_idx] - coords[water_switch]));
            }
            return mi;
        }

        double mi = std::numeric_limits<double>::infinity();

        for (int atom1 : groups[group1].atoms) {
            const int atom1_idx = atom1 - 1;

            for (int atom2 : groups[group2].atoms) {
                const int atom2_idx = atom2 - 1;
                mi = std::min(mi, norm2(coords[atom1_idx] - coords[atom2_idx]));
            }
        }
        return mi;
    };

    auto group_is_active = [&](int group) {
        const int switch_atom = groups[group].iswitch - 1;
        return switch_atom >= 0 && switch_atom < ctx.n_atoms && !excluded[switch_atom];
    };

    const double lrf_cutoff2 = ctx.md.lrf_cutoff * ctx.md.lrf_cutoff;

    exact_calculation_groups_.clear();
    lrf_calculation_groups_.clear();

    for (int i = 0; i < n_groups; i++) {
        if (!group_is_active(i)) continue;
        for (int j = i; j < n_groups; j++) {
            if (!group_is_active(j)) continue;

            const double distance2 = group_distance2(i, j);
            const double cutoff = normal_cutoff(i, j);
            const double cutoff2 = cutoff * cutoff;

            if (distance2 <= cutoff2) {
                // need to calculate each pair
                exact_calculation_groups_.push_back({i, j});
            } else if (distance2 <= lrf_cutoff2) {
                // need to use lrf
                lrf_calculation_groups_.push_back({i, j});
            } else {
                // ignore
                continue;
            }
        }
    }
    init_exact_atom_pairs(ctx);

    printf(
        "[Non Bonded Force] there are %d exact atom pairs. Total number of atoms: %d. Exact atom pairs ratio is %f\n",
        (int)exact_atom_pairs_.size(), ctx.n_atoms, 1.0 * exact_atom_pairs_.size() / ctx.n_atoms / ctx.n_atoms);
}

void CpuNonbondedForce::init_exact_atom_pairs(Context& ctx) {
    const auto& groups = ctx.charge_group_config.charge_groups;
    const int n_groups = groups.size();

    const auto* atom_to_group = data_.atom_to_group->cpu_data_p;

    // Dense lookup table for exact charge-group pairs.
    std::vector<uint8_t> is_exact_group_pair(static_cast<size_t>(n_groups) * n_groups, 0);

    for (const auto& pair : exact_calculation_groups_) {
        is_exact_group_pair[static_cast<size_t>(pair.first) * n_groups + pair.second] = 1;

        is_exact_group_pair[static_cast<size_t>(pair.second) * n_groups + pair.first] = 1;
    }

    exact_atom_pairs_.clear();

    const int* atom_indices = data_.atom_idx->cpu_data_p;
    const uint8_t* categories = data_.category->cpu_data_p;
    const int* q_states = data_.q_state->cpu_data_p;

    constexpr uint8_t Q = static_cast<uint8_t>(AtomCategory::Q);

    const double rcq2 = ctx.md.q_atom * ctx.md.q_atom;
    const coord_t& solute_center = ctx.topo.solute_center;
    const auto& config = ctx.charge_group_config;
    const coord_t* coords = ctx.coords->cpu_data_p;

    auto inside_rcq = [&](int atom, uint8_t category) {
        if (category == static_cast<uint8_t>(AtomCategory::W)) {
            const int water_offset = atom - ctx.n_atoms_solute;
            const int water_oxygen = ctx.n_atoms_solute + (water_offset / 3) * 3;
            return norm2(coords[water_oxygen] - solute_center) <= rcq2;
        } else if (category == static_cast<uint8_t>(AtomCategory::P)) {
            const int group = atom_to_group[atom];
            if (group < 0) return false;
            if (config.iuse_switch_atom == 1) {
                const int switch_atom = groups[group].iswitch - 1;
                return norm2(coords[switch_atom] - solute_center) <= rcq2;
            }
            for (int atom_1based : groups[group].atoms) {
                const int group_atom = atom_1based - 1;
                if (norm2(coords[group_atom] - solute_center) <= rcq2) {
                    return true;
                }
            }
            return false;
        }
        return true;
    };

    for (int slot1 = 0; slot1 < data_.n_total; slot1++) {
        const int atom1 = atom_indices[slot1];
        if (atom1 < 0) continue;

        const bool atom1_is_q = categories[slot1] == Q;
        for (int slot2 = slot1 + 1; slot2 < data_.n_total; slot2++) {
            const int atom2 = atom_indices[slot2];
            if (atom2 < 0 || atom1 == atom2) continue;

            const bool atom2_is_q = categories[slot2] == Q;

            bool calculate_directly = false;

            if (atom1_is_q && atom2_is_q) {
                calculate_directly = true;
            } else if (atom1_is_q || atom2_is_q) {
                /*
                 * LRF excludes Q atoms. Preserve the current QGPU behavior by
                 * calculating every Q-containing interaction directly.
                 */
                const int environment_atom = atom1_is_q ? atom2 : atom1;
                const uint8_t environment_category = atom1_is_q ? categories[slot2] : categories[slot1];
                calculate_directly = inside_rcq(environment_atom, environment_category);

            } else {
                const int group1 = atom_to_group[atom1];
                const int group2 = atom_to_group[atom2];

                if (group1 < 0 || group2 < 0) {
                    continue;
                }

                calculate_directly = is_exact_group_pair[static_cast<size_t>(group1) * n_groups + group2] != 0;
            }

            if (!calculate_directly) continue;

            exact_atom_pairs_.push_back({slot1, slot2});
        }
    }
}

void CpuNonbondedForce::calc_direct_pair(Context& ctx, int slot1, int slot2) {
    const auto& atom_idxs = data_.atom_idx->cpu_data_p;
    const auto& coords = ctx.coords->cpu_data_p;
    auto& dvelocities = ctx.dvelocities->cpu_data_p;
    int sz = data_.n_total;

    const int atom1 = atom_idxs[slot1];
    const int atom2 = atom_idxs[slot2];

    if (atom1 == -1 || atom2 == -1) return;

    const auto& atom1_type = data_.category->cpu_data_p[slot1];
    const int atom1_state = data_.q_state->cpu_data_p[slot1];
    const real_t atom1_charge = data_.atom_charge->cpu_data_p[slot1];
    const vdw_atom_param_t& atom1_vdw = data_.atom_vdw->cpu_data_p[slot1];

    const auto& atom2_type = data_.category->cpu_data_p[slot2];
    const int atom2_state = data_.q_state->cpu_data_p[slot2];
    const auto& bond_type = get_bond_type(ctx.n_atoms_solute, ctx.LJ_matrix->cpu_data_p, atom1, atom1_type, atom2, atom2_type);
    const real_t atom2_charge = data_.atom_charge->cpu_data_p[slot2];
    const vdw_atom_param_t atom2_vdw = data_.atom_vdw->cpu_data_p[slot2];

    if (bond_type == BondType::Bond23) return;
    if (atom1_type == static_cast<uint8_t>(AtomCategory::Q) && atom2_type == static_cast<uint8_t>(AtomCategory::Q) && atom1_state != atom2_state) {
        return;
    }

    real_t dx = coords[atom2].x - coords[atom1].x;
    real_t dy = coords[atom2].y - coords[atom1].y;
    real_t dz = coords[atom2].z - coords[atom1].z;
    real_t dis2 = dx * dx + dy * dy + dz * dz;
    real_t inv_dis2 = static_cast<real_t>(1.0) / dis2;
    real_t inv_dis = sqrt(inv_dis2);

    real_t qij = atom1_charge * atom2_charge;
    bool is_14 = (bond_type == BondType::Bond14);
    real_t scaling = is_14 ? ctx.topo.el14_scale : 1;
    real_t2 pair = is_14 ? combine_vdw(ctx.topo.vdw_rule, atom1_vdw.aii_14, atom1_vdw.bii_14, atom2_vdw.aii_14, atom2_vdw.bii_14) : combine_vdw(ctx.topo.vdw_rule, atom1_vdw.aii_normal, atom1_vdw.bii_normal, atom2_vdw.aii_normal, atom2_vdw.bii_normal);

    auto [vel, dvel] = calc_electrostatic(qij * scaling, ctx.topo.coulomb_constant, inv_dis);
    auto [vvdw, dvvdw] = calc_vdw(pair, inv_dis);

    real_t lambda = std::min(data_.atom_lambdas->cpu_data_p[slot1], data_.atom_lambdas->cpu_data_p[slot2]);

    real_t dva = (dvel + dvvdw) * inv_dis * lambda;

    add_force(dvelocities[atom1].x, -dva * dx);
    add_force(dvelocities[atom1].y, -dva * dy);
    add_force(dvelocities[atom1].z, -dva * dz);

    add_force(dvelocities[atom2].x, dva * dx);
    add_force(dvelocities[atom2].y, dva * dy);
    add_force(dvelocities[atom2].z, dva * dz);

    // Accumulate energy
    accumulate_energy(ctx, vel, vvdw, atom1_type, atom2_type, atom1_state, atom2_state);
}

void CpuNonbondedForce::calc_exact_pairs(Context& ctx) {
    for (const auto& pair : exact_atom_pairs_) {
        calc_direct_pair(ctx, pair.first, pair.second);
    }
}

void CpuNonbondedForce::init_backend(Context& ctx) {
    non_q_slot_by_atom_.assign(ctx.n_atoms, -1);

    const int* atom_indices = data_.atom_idx->cpu_data_p;
    const uint8_t* categories = data_.category->cpu_data_p;

    constexpr uint8_t P = static_cast<uint8_t>(AtomCategory::P);
    constexpr uint8_t W = static_cast<uint8_t>(AtomCategory::W);

    for (int slot = 0; slot < data_.n_total; ++slot) {
        const int atom = atom_indices[slot];
        if (atom < 0) continue;

        if (categories[slot] == P || categories[slot] == W) {
            non_q_slot_by_atom_[atom] = slot;
        }
    }
}

void CpuNonbondedForce::init_lrf_coefficients(Context& ctx) {
    const auto& groups = ctx.charge_group_config.charge_groups;

    const int n_groups = static_cast<int>(groups.size());

    const coord_t* coords = ctx.coords->cpu_data_p;

    lrf_coefficients_.assign(n_groups, LrfCoefficients{});

    for (int group = 0; group < n_groups; group++) {
        coord_t center = {};
        const auto& group_atoms = groups[group].atoms;
        for (int atom_1based : group_atoms) {
            const int atom = atom_1based - 1;

            center = center + coords[atom];
        }
        const double inverse_count = 1.0 / group_atoms.size();
        center = center * inverse_count;

        lrf_coefficients_[group].center = center;
    }

    const real_t* charges = data_.atom_charge->cpu_data_p;

    auto accumulate_group_into_target = [&](int source_group, int target_group) {
        LrfCoefficients& target = lrf_coefficients_[target_group];

        for (int atom_1based : groups[source_group].atoms) {
            const int atom = atom_1based - 1;
            const int slot = non_q_slot_by_atom_[atom];
            if (slot < 0) continue;
            accumulate_lrf_source(target, coords[atom], charges[slot]);
        }
    };

    for (const auto& pair : lrf_calculation_groups_) {
        const int group1 = pair.first;
        const int group2 = pair.second;
        accumulate_group_into_target(group1, group2);
        accumulate_group_into_target(group2, group1);
    }
}

void CpuNonbondedForce::calc_lrf(Context& ctx) {
    const int* atom_indices = data_.atom_idx->cpu_data_p;
    const uint8_t* categories = data_.category->cpu_data_p;
    const real_t* charges = data_.atom_charge->cpu_data_p;
    auto* atom_to_group = data_.atom_to_group->cpu_data_p;

    const coord_t* coords = ctx.coords->cpu_data_p;
    auto* dvelocities = ctx.dvelocities->cpu_data_p;

    constexpr uint8_t P = static_cast<uint8_t>(AtomCategory::P);
    constexpr uint8_t W = static_cast<uint8_t>(AtomCategory::W);

    const double coulomb_constant = ctx.topo.coulomb_constant;
    double lrf_energy = 0.0;

    for (int slot = 0; slot < data_.n_total; slot++) {
        const int atom = atom_indices[slot];
        if (atom < 0) continue;

        // LRF is applied only to non-Q atoms
        if (categories[slot] != P && categories[slot] != W) {
            continue;
        }

        const int group = atom_to_group[atom];
        if (group < 0) continue;

        const LrfCoefficients& lrf = lrf_coefficients_[group];

        coord_t d = lrf.center - coords[atom];
        double d_array[3] = {d.x, d.y, d.z};

        double potential = lrf.phi0;

        for (int a = 0; a < 3; a++) {
            potential += lrf.phi1[a] * d_array[a];
        }

        for (int a = 0; a < 3; a++) {
            for (int b = 0; b < 3; b++) {
                potential += 0.5 * lrf.phi2[a * 3 + b] * d_array[a] * d_array[b];
            }
        }

        double df[3] = {lrf.phi1[0], lrf.phi1[1], lrf.phi1[2]};
        for (int a = 0; a < 3; a++) {
            for (int b = 0; b < 3; b++) {
                df[a] += lrf.phi2[a * 3 + b] * d_array[b];
            }
        }

        for (int a = 0; a < 3; a++) {
            for (int b = 0; b < 3; b++) {
                for (int c = 0; c < 3; c++) {
                    df[a] += 0.5 * lrf.phi3[(a * 3 + b) * 3 + c] * d_array[b] * d_array[c];
                }
            }
        }

        const double charge = static_cast<double>(charges[slot]);

        lrf_energy += 0.5 * coulomb_constant * charge * potential;

        const double gradient_scale = -coulomb_constant * charge;

        add_force(dvelocities[atom].x, gradient_scale * df[0]);
        add_force(dvelocities[atom].y, gradient_scale * df[1]);
        add_force(dvelocities[atom].z, gradient_scale * df[2]);
    }

    add_energy(ctx.energy.host()[E_LRF], lrf_energy);
}

void CpuNonbondedForce::calc(Context& ctx) {
    if (!ctx.md.lrf || ctx.md.non_bond == 0) {
        calc_all_direct_pairs(ctx);
        return;
    }

    // 1. Examine each pair of charge groups
    if (ctx.step == ctx.md.steps || ctx.step % ctx.md.non_bond == 0) {
        init_calculation_groups(ctx);
        init_lrf_coefficients(ctx);
    }
    calc_exact_pairs(ctx);
    calc_lrf(ctx);
}