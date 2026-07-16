#include "bonded_force.h"

namespace {
inline double deg2rad(double d) {
    return d * (M_PI / 180.0);
}

}  // namespace

void BondedForce::init(Context& ctx, const ParseResult& parsed, const ShakeData& shake_data) {
    build_bonds(ctx, parsed, shake_data);
    build_angles(ctx, parsed, shake_data);
    build_torsions(ctx, parsed, shake_data);
    build_impropers(ctx, parsed, shake_data);
    init_backend(ctx);
}

void BondedForce::build_bonds(Context& ctx, const ParseResult& parsed, const ShakeData& shake_data) {
    const bool run_gpu = ctx.command_info.requested_gpu;

    std::vector<bond_idx_t> ids;
    std::vector<dparam2_t> params;
    std::vector<int> eslot;

    for (int i = 0; i < parsed.bonds.size(); i++) {
        const bond_t& bond = parsed.bonds[i];
        if (shake_data.contains(bond.ai, bond.aj)) {
            continue;
        }

        const cbond_t& c = parsed.cbonds[bond.code - 1];
        ids.push_back({bond.ai - 1, bond.aj - 1});
        params.push_back({c.kb, c.b0});  // {kb, b0}
        eslot.push_back(i < parsed.n_bonds_solute ? E_BOND_P_BOND : E_BOND_W_BOND);
    }

    data_.bond.n = static_cast<int>(ids.size());
    data_.bond.ids = HostDeviceBuffer<bond_idx_t>::from_vector(ids, run_gpu);
    data_.bond.params = HostDeviceBuffer<dparam2_t>::from_vector(params, run_gpu);
    data_.bond.eslot = HostDeviceBuffer<int>::from_vector(eslot, run_gpu);
}
void BondedForce::build_angles(Context& ctx, const ParseResult& parsed, const ShakeData& shake_data) {
    const bool run_gpu = ctx.command_info.requested_gpu;

    std::vector<angle_idx_t> ids;
    std::vector<dparam2_t> params;
    std::vector<int> eslot;

    for (int i = 0; i < parsed.angles.size(); i++) {
        const angle_t& angle = parsed.angles[i];
        if (shake_data.contains(angle.ai, angle.ak)) {
            continue;
        }

        const cangle_t& c = parsed.cangles[angle.code - 1];
        ids.push_back({angle.ai - 1, angle.aj - 1, angle.ak - 1});
        params.push_back({c.kth, deg2rad(c.th0)});  // {kth, th0(rad)}
        eslot.push_back(i < parsed.n_angles_solute ? E_BOND_P_ANGLE : E_BOND_W_ANGLE);
    }

    data_.angle.n = static_cast<int>(ids.size());
    data_.angle.ids = HostDeviceBuffer<angle_idx_t>::from_vector(ids, run_gpu);
    data_.angle.params = HostDeviceBuffer<dparam2_t>::from_vector(params, run_gpu);
    data_.angle.eslot = HostDeviceBuffer<int>::from_vector(eslot, run_gpu);
}
void BondedForce::build_torsions(Context& ctx, const ParseResult& parsed, const ShakeData& shake_data) {
    const bool run_gpu = ctx.command_info.requested_gpu;

    std::vector<dihe_idx_t> ids;
    std::vector<torsion_param_t> params;
    std::vector<int> eslot;

    for (int i = 0; i < parsed.torsions.size(); i++) {
        
        const torsion_t& torsion = parsed.torsions[i];

        const ctorsion_t& c = parsed.ctorsions[torsion.code - 1];
        ids.push_back({torsion.ai - 1, torsion.aj - 1,
                       torsion.ak - 1, torsion.al - 1});
        params.push_back({static_cast<real_t>(c.k),
                          static_cast<real_t>(c.n),
                          static_cast<real_t>(deg2rad(c.d)),
                          static_cast<real_t>(c.paths)});  // {k, n, d(rad), paths}
        eslot.push_back(i < parsed.n_torsions_solute ? E_BOND_P_TOR : E_BOND_W_TOR);
    }

    data_.torsion.n = static_cast<int>(ids.size());
    data_.torsion.ids = HostDeviceBuffer<dihe_idx_t>::from_vector(ids, run_gpu);
    data_.torsion.params = HostDeviceBuffer<torsion_param_t>::from_vector(params, run_gpu);
    data_.torsion.eslot = HostDeviceBuffer<int>::from_vector(eslot, run_gpu);
}
void BondedForce::build_impropers(Context& ctx, const ParseResult& parsed, const ShakeData& shake_data) {
    const bool run_gpu = ctx.command_info.requested_gpu;

    std::vector<dihe_idx_t> ids;
    std::vector<dparam2_t> params;
    std::vector<int> eslot;

    for (int i = 0; i < parsed.impropers.size(); i++) {

        const improper_t& improper = parsed.impropers[i];

        const cimproper_t& c = parsed.cimpropers[improper.code - 1];
        ids.push_back({improper.ai - 1, improper.aj - 1,
                       improper.ak - 1, improper.al - 1});
        params.push_back({c.k, deg2rad(c.phi0)});  // {k, phi0(rad)}
        eslot.push_back(i < parsed.n_impropers_solute ? E_BOND_P_IMP : E_BOND_W_IMP);
    }

    data_.improper.n = static_cast<int>(ids.size());
    data_.improper.ids = HostDeviceBuffer<dihe_idx_t>::from_vector(ids, run_gpu);
    data_.improper.params = HostDeviceBuffer<dparam2_t>::from_vector(params, run_gpu);
    data_.improper.eslot = HostDeviceBuffer<int>::from_vector(eslot, run_gpu);
}
