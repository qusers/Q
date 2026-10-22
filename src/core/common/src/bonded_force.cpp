#include "bonded_force.h"

namespace {
inline double deg2rad(double d) {
    return d * (M_PI / 180.0);
}

using BondKey = std::pair<int, int>;
using AngleKey = std::array<int, 3>;
using DihedralKey = std::array<int, 4>;

BondKey bond_key(int i, int j) {
    return {std::min(i, j), std::max(i, j)};
}

AngleKey angle_key(int i, int j, int k) {
    return {std::min(i, k), j, std::max(i, k)};
}

DihedralKey torsion_key(int i, int j, int k, int l) {
    const DihedralKey forward{i, j, k, l};
    const DihedralKey reverse{l, k, j, i};
    return std::min(forward, reverse);
}

BondKey improper_key(int j, int k) {
    return bond_key(j, k);
}

size_t terms_per_state(size_t total, int nstates) {
    if (total == 0) return 0;

    if (nstates <= 0 || total % static_cast<size_t>(nstates) != 0) {
        throw std::runtime_error("Invalid state-dependent bonded interaction layout.");
    }

    return total / static_cast<size_t>(nstates);
}

int q_energy_slot(int state, EQSlot slot) {
    return EnergyBuffer::eq_index(ENERGY_FIXED_COUNT, state, slot);
}

}  // namespace

void BondedForce::init(Context& ctx, const ParseResult& parsed, const ConstraintData& shake_data, const std::vector<special_pair_term_t>& special_pairs) {
    if (!parsed.q_shakes.empty()) {
        throw std::runtime_error("State-dependent Q constraints are not implemented.");
    }

    for (const auto& bond : parsed.q_bonds) {
        if (shake_data.contains(bond.ai, bond.aj)) {
            throw std::runtime_error("A redefined Q-bond is constrained.");
        }
    }

    for (const auto& angle : parsed.q_angles) {
        if (shake_data.contains(angle.ai, angle.ak)) {
            throw std::runtime_error("A redefined Q-angle has constrained outer atoms.");
        }
    }

    build_bonds(ctx, parsed, shake_data);
    build_angles(ctx, parsed, shake_data);
    build_torsions(ctx, parsed, shake_data);
    build_impropers(ctx, parsed, shake_data);
    build_special_pairs(ctx, special_pairs);
    init_backend(ctx);
}

void BondedForce::build_bonds(Context& ctx, const ParseResult& parsed, const ConstraintData& shake_data) {
    const bool run_gpu = ctx.command_info.requested_gpu;

    std::vector<bond_idx_t> ids;
    std::vector<cbond_t> params;
    std::vector<int> eslot;
    std::vector<int> state;

    // Include code-zero Q records: they also replace topology terms
    std::set<BondKey> replaced;
    for (const auto& bond : parsed.q_bonds) {
        replaced.insert(bond_key(bond.ai, bond.aj));
    }

    // Normal bonds
    for (int i = 0; i < parsed.bonds.size(); i++) {
        const bond_t& bond = parsed.bonds[i];
        if (bond.code == 0 || replaced.count(bond_key(bond.ai, bond.aj)) || shake_data.contains(bond.ai, bond.aj)) {
            continue;
        }

        const cbond_t& c = parsed.cbonds[bond.code - 1];
        ids.push_back({bond.ai - 1, bond.aj - 1});
        params.push_back(c);
        eslot.push_back(i < parsed.n_bonds_solute ? E_BOND_P_BOND : E_BOND_W_BOND);
        state.push_back(-1);
    }

    // Q bonds
    int nstates = ctx.n_lambdas();
    const size_t nq = terms_per_state(parsed.q_bonds.size(), nstates);
    for (int i = 0; i < parsed.q_bonds.size(); i++) {
        const auto& bond = parsed.q_bonds[i];
        if (bond.code == 0) {
            continue;
        }

        int s = i / nq;

        const auto& c = parsed.q_cbonds.at(bond.code);
        ids.push_back({bond.ai - 1, bond.aj - 1});
        params.push_back(c);
        state.push_back(s);
        eslot.push_back(q_energy_slot(s, EQ_BOND_BOND));
    }

    data_.bond.n = static_cast<int>(ids.size());
    data_.bond.ids = HostDeviceBuffer<bond_idx_t>::from_vector(ids, run_gpu);
    data_.bond.params = HostDeviceBuffer<cbond_t>::from_vector(params, run_gpu);
    data_.bond.eslot = HostDeviceBuffer<int>::from_vector(eslot, run_gpu);
    data_.bond.state = HostDeviceBuffer<int>::from_vector(state, run_gpu);
}
void BondedForce::build_angles(Context& ctx, const ParseResult& parsed, const ConstraintData& shake_data) {
    const bool run_gpu = ctx.command_info.requested_gpu;

    std::vector<angle_idx_t> ids;
    std::vector<dparam2_t> params;
    std::vector<int> eslot;
    std::vector<int> state;
    std::vector<morse_coupling_t> couplings;
    const int nstates = ctx.n_lambdas();
    const size_t nq = terms_per_state(parsed.q_angles.size(), nstates);
    const size_t nbonds = terms_per_state(parsed.q_bonds.size(), nstates);

    // Original Q-angle row -> original Q-bond row.
    // Both indices are 0-based here; -1 means uncoupled.
    std::vector<int> bond_for_angle(nq, -1);

    for (const auto& couple : parsed.q_angcouples) {
        bond_for_angle.at(couple.angle_index - 1) =
            couple.bond_index - 1;
    }

    std::set<AngleKey> replaced;
    for (const auto& angle : parsed.q_angles) {
        replaced.insert(angle_key(angle.ai, angle.aj, angle.ak));
    }

    // normal angle
    for (int i = 0; i < parsed.angles.size(); i++) {
        const angle_t& angle = parsed.angles[i];
        if (angle.code == 0 || replaced.count(angle_key(angle.ai, angle.aj, angle.ak)) ||
            shake_data.contains(angle.ai, angle.ak)) {
            continue;
        }

        const cangle_t& c = parsed.cangles[angle.code - 1];
        ids.push_back({angle.ai - 1, angle.aj - 1, angle.ak - 1});
        params.push_back({c.kth, deg2rad(c.th0)});  // {kth, th0(rad)}
        eslot.push_back(i < parsed.n_angles_solute ? E_BOND_P_ANGLE : E_BOND_W_ANGLE);
        state.push_back(-1);
        couplings.emplace_back();  // i = j = -1: uncoupled
    }

    for (int i = 0; i < parsed.q_angles.size(); i++) {
        const auto& angle = parsed.q_angles[i];
        if (angle.code == 0) {
            continue;
        }

        int s = i / nq;
        const size_t angle_index = i % nq;

        const auto& c = parsed.q_cangles.at(angle.code);

        morse_coupling_t coupling{};
        const int bond_index = bond_for_angle.at(angle_index);
        if (bond_index >= 0) {
            const auto& bond = parsed.q_bonds.at(static_cast<size_t>(s) * nbonds + bond_index);

            if (bond.code <= 0) {
                throw std::runtime_error("A coupled Q-angle references an absent bond.");
            }

            const auto& bond_param = parsed.q_cbonds.at(bond.code);
            if (!bond_param.is_morse) {
                throw std::runtime_error("A coupled Q-angle must reference a Morse bond.");
            }

            coupling.i = bond.ai - 1;
            coupling.j = bond.aj - 1;
            coupling.a = bond_param.a;
            coupling.r0 = bond_param.b0;
        }

        ids.push_back({angle.ai - 1, angle.aj - 1, angle.ak - 1});
        params.push_back({c.kth, deg2rad(c.th0)});
        state.push_back(s);
        eslot.push_back(q_energy_slot(s, EQ_BOND_ANGLE));
        couplings.push_back(coupling);
    }

    data_.angle.n = static_cast<int>(ids.size());
    data_.angle.ids = HostDeviceBuffer<angle_idx_t>::from_vector(ids, run_gpu);
    data_.angle.params = HostDeviceBuffer<dparam2_t>::from_vector(params, run_gpu);
    data_.angle.eslot = HostDeviceBuffer<int>::from_vector(eslot, run_gpu);
    data_.angle.state = HostDeviceBuffer<int>::from_vector(state, run_gpu);
    data_.angle.couplings = HostDeviceBuffer<morse_coupling_t>::from_vector(couplings, run_gpu);
}
void BondedForce::build_torsions(Context& ctx, const ParseResult& parsed, const ConstraintData& shake_data) {
    const bool run_gpu = ctx.command_info.requested_gpu;

    std::vector<dihe_idx_t> ids;
    std::vector<torsion_param_t> params;
    std::vector<int> eslot;
    std::vector<int> state;
    std::vector<morse_coupling_t> couplings;

    const int nstates = ctx.n_lambdas();
    const size_t nq = terms_per_state(parsed.q_torsions.size(), nstates);
    const size_t nbonds = terms_per_state(parsed.q_bonds.size(), nstates);

    // Original Q-torsion row -> original Q-bond row.
    // Both indices are 0-based here; -1 means uncoupled.
    std::vector<int> bond_for_torsion(nq, -1);

    for (const auto& couple : parsed.q_torcouples) {
        bond_for_torsion.at(couple.torsion_index - 1) = couple.bond_index - 1;
    }

    std::set<DihedralKey> replaced;
    for (const auto& torsion : parsed.q_torsions) {
        replaced.insert(torsion_key(torsion.ai, torsion.aj, torsion.ak, torsion.al));
    }

    for (int i = 0; i < parsed.torsions.size(); i++) {
        const torsion_t& torsion = parsed.torsions[i];

        if (torsion.code == 0 || replaced.count(torsion_key(torsion.ai, torsion.aj, torsion.ak, torsion.al))) {
            continue;
        }

        const ctorsion_t& c = parsed.ctorsions[torsion.code - 1];
        ids.push_back({torsion.ai - 1, torsion.aj - 1,
                       torsion.ak - 1, torsion.al - 1});
        params.push_back({static_cast<real_t>(c.k),
                          static_cast<real_t>(c.n),
                          static_cast<real_t>(deg2rad(c.d)),
                          static_cast<real_t>(c.paths)});  // {k, n, d(rad), paths}
        eslot.push_back(i < parsed.n_torsions_solute ? E_BOND_P_TOR : E_BOND_W_TOR);
        state.push_back(-1);
        couplings.emplace_back();
    }

    // Q torsion
    for (int i = 0; i < parsed.q_torsions.size(); i++) {
        const auto& torsion = parsed.q_torsions[i];
        if (torsion.code == 0) {
            continue;
        }

        int s = i / nq;
        const size_t torsion_index = i % nq;

        const auto& c = parsed.q_ctorsions.at(torsion.code);

        morse_coupling_t coupling{};
        const int bond_index = bond_for_torsion.at(torsion_index);

        if (bond_index >= 0) {
            const auto& bond = parsed.q_bonds.at(static_cast<size_t>(s) * nbonds + bond_index);

            if (bond.code <= 0) {
                throw std::runtime_error("A coupled Q-torsion references an absent bond.");
            }

            const auto& bond_param = parsed.q_cbonds.at(bond.code);
            if (!bond_param.is_morse) {
                throw std::runtime_error("A coupled Q-torsion must reference a Morse bond.");
            }

            coupling.i = bond.ai - 1;
            coupling.j = bond.aj - 1;
            coupling.a = bond_param.a;
            coupling.r0 = bond_param.b0;
        }

        ids.push_back({torsion.ai - 1, torsion.aj - 1, torsion.ak - 1, torsion.al - 1});
        params.push_back({static_cast<real_t>(c.k), static_cast<real_t>(c.n), static_cast<real_t>(deg2rad(c.d)), static_cast<real_t>(1.0)});
        state.push_back(s);
        eslot.push_back(q_energy_slot(s, EQ_BOND_TOR));
        couplings.push_back(coupling);
    }

    data_.torsion.n = static_cast<int>(ids.size());
    data_.torsion.ids = HostDeviceBuffer<dihe_idx_t>::from_vector(ids, run_gpu);
    data_.torsion.params = HostDeviceBuffer<torsion_param_t>::from_vector(params, run_gpu);
    data_.torsion.eslot = HostDeviceBuffer<int>::from_vector(eslot, run_gpu);
    data_.torsion.state = HostDeviceBuffer<int>::from_vector(state, run_gpu);
    data_.torsion.couplings = HostDeviceBuffer<morse_coupling_t>::from_vector(couplings, run_gpu);
}
void BondedForce::build_impropers(Context& ctx, const ParseResult& parsed, const ConstraintData& shake_data) {
    const bool run_gpu = ctx.command_info.requested_gpu;

    std::vector<dihe_idx_t> ids;
    std::vector<dparam2_t> params;
    std::vector<int> eslot;
    std::vector<int> state;
    std::vector<morse_coupling_t> couplings;

    const int nstates = ctx.n_lambdas();
    const size_t nq = terms_per_state(parsed.q_impropers.size(), nstates);
    const size_t nbonds = terms_per_state(parsed.q_bonds.size(), nstates);

    // Original Q-improper row -> original Q-bond row.
    // Both indices are 0-based here.
    std::vector<int> bond_for_improper(nq, -1);
    std::vector<bool> inverted_for_improper(nq, false);

    for (const auto& couple : parsed.q_imprcouples) {
        bond_for_improper.at(couple.improper_index - 1) = couple.bond_index - 1;
        inverted_for_improper.at(couple.improper_index - 1) = couple.inverted;
    }

    /*
    todo: The replacement rule is different from the force field.
        Now support the default replacement rule.
        If want to support all, should read FF_TYPE in top file, to get the force field version.
    */
    std::set<BondKey> replaced;
    for (const auto& improper : parsed.q_impropers) {
        replaced.insert(improper_key(improper.aj, improper.ak));
    }

    for (int i = 0; i < parsed.impropers.size(); i++) {
        const improper_t& improper = parsed.impropers[i];

        if (improper.code == 0 || replaced.count(improper_key(improper.aj, improper.ak))) {
            continue;
        }

        const cimproper_t& c = parsed.cimpropers[improper.code - 1];
        ids.push_back({improper.ai - 1, improper.aj - 1,
                       improper.ak - 1, improper.al - 1});
        params.push_back({c.k, deg2rad(c.phi0)});  // {k, phi0(rad)}
        eslot.push_back(i < parsed.n_impropers_solute ? E_BOND_P_IMP : E_BOND_W_IMP);
        state.push_back(-1);
        couplings.emplace_back();
    }

    for (int i = 0; i < parsed.q_impropers.size(); i++) {
        const auto& improper = parsed.q_impropers[i];

        if (improper.code == 0) {
            continue;
        }

        int s = i / nq;
        const size_t improper_index = i % nq;

        const auto& c = parsed.q_cimpropers.at(improper.code);

        morse_coupling_t coupling{};
        const int bond_index = bond_for_improper.at(improper_index);

        if (bond_index >= 0) {
            const auto& bond = parsed.q_bonds.at(static_cast<size_t>(s) * nbonds + bond_index);

            if (bond.code <= 0) {
                throw std::runtime_error("A coupled Q-improper references an absent bond.");
            }

            const auto& bond_param = parsed.q_cbonds.at(bond.code);
            if (!bond_param.is_morse) {
                throw std::runtime_error("A coupled Q-improper must reference a Morse bond.");
            }

            coupling.i = bond.ai - 1;
            coupling.j = bond.aj - 1;
            coupling.a = bond_param.a;
            coupling.r0 = bond_param.b0;
            coupling.inverted = inverted_for_improper.at(improper_index);
        }

        ids.push_back({improper.ai - 1, improper.aj - 1, improper.ak - 1, improper.al - 1});
        params.push_back({c.k, deg2rad(c.phi0)});
        state.push_back(s);
        eslot.push_back(q_energy_slot(s, EQ_BOND_IMP));
        couplings.push_back(coupling);
    }

    data_.improper.n = static_cast<int>(ids.size());
    data_.improper.ids = HostDeviceBuffer<dihe_idx_t>::from_vector(ids, run_gpu);
    data_.improper.params = HostDeviceBuffer<dparam2_t>::from_vector(params, run_gpu);
    data_.improper.eslot = HostDeviceBuffer<int>::from_vector(eslot, run_gpu);
    data_.improper.state = HostDeviceBuffer<int>::from_vector(state, run_gpu);
    data_.improper.couplings = HostDeviceBuffer<morse_coupling_t>::from_vector(couplings, run_gpu);
}

void BondedForce::build_special_pairs(Context& ctx, const std::vector<special_pair_term_t>& special_pairs) {
    data_.special_pairs.n = static_cast<int>(special_pairs.size());

    data_.special_pairs.terms = HostDeviceBuffer<special_pair_term_t>::from_vector(special_pairs, ctx.command_info.requested_gpu);
}