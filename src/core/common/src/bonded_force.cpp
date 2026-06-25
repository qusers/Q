#pragma once
#include "bonded_force.h"

namespace {
inline double deg2rad(double d) {
    return d * (M_PI / 180.0);
}

}  // namespace

void BondedForce::init(Context& ctx) {
    build_bonds(ctx);
    build_angles(ctx);
    build_torsions(ctx);
    build_impropers(ctx);
    init_backend(ctx);
}

void BondedForce::build_bonds(Context& ctx) {
    const bool run_gpu = ctx.command_info.requested_gpu;
    auto* bonds = ctx.bonds->cpu_data_p;
    auto* cbonds = ctx.cbonds->cpu_data_p;

    std::vector<bond_idx_t> ids;
    std::vector<dparam2_t> params;
    std::vector<int> eslot;

    ids.reserve(ctx.n_bonds);
    params.reserve(ctx.n_bonds);
    eslot.reserve(ctx.n_bonds);

    for (int i = 0; i < ctx.n_bonds; i++) {
        const cbond_t& c = cbonds[bonds[i].code - 1];
        ids.push_back({bonds[i].ai - 1, bonds[i].aj - 1});
        params.push_back({c.kb, c.b0});  // {kb, b0}
        eslot.push_back(i < ctx.n_bonds_solute ? bonded_p_slot() : bonded_w_slot());
    }

    data_.bond.n = static_cast<int>(ids.size());
    data_.bond.ids = HostDeviceBuffer<bond_idx_t>::from_vector(ids, run_gpu);
    data_.bond.params = HostDeviceBuffer<dparam2_t>::from_vector(params, run_gpu);
    data_.bond.eslot = HostDeviceBuffer<int>::from_vector(eslot, run_gpu);
}
void BondedForce::build_angles(Context& ctx) {
    const bool run_gpu = ctx.command_info.requested_gpu;
    auto* angles = ctx.angles->cpu_data_p;
    auto* cangles = ctx.cangles->cpu_data_p;

    std::vector<angle_idx_t> ids;
    std::vector<dparam2_t> params;
    std::vector<int> eslot;
    ids.reserve(ctx.n_angles);
    params.reserve(ctx.n_angles);
    eslot.reserve(ctx.n_angles);

    for (int i = 0; i < ctx.n_angles; i++) {
        const cangle_t& c = cangles[angles[i].code - 1];
        ids.push_back({angles[i].ai - 1, angles[i].aj - 1, angles[i].ak - 1});
        params.push_back({c.kth, deg2rad(c.th0)});  // {kth, th0(rad)}
        eslot.push_back(i < ctx.n_angles_solute ? bonded_p_slot() : bonded_w_slot());
    }

    data_.angle.n = static_cast<int>(ids.size());
    data_.angle.ids = HostDeviceBuffer<angle_idx_t>::from_vector(ids, run_gpu);
    data_.angle.params = HostDeviceBuffer<dparam2_t>::from_vector(params, run_gpu);
    data_.angle.eslot = HostDeviceBuffer<int>::from_vector(eslot, run_gpu);
}
void BondedForce::build_torsions(Context& ctx) {
    const bool run_gpu = ctx.command_info.requested_gpu;
    auto* torsions = ctx.torsions->cpu_data_p;
    auto* ctorsions = ctx.ctorsions->cpu_data_p;

    std::vector<dihe_idx_t> ids;
    std::vector<torsion_param_t> params;
    std::vector<int> eslot;
    ids.reserve(ctx.n_torsions);
    params.reserve(ctx.n_torsions);
    eslot.reserve(ctx.n_torsions);

    for (int i = 0; i < ctx.n_torsions; i++) {
        const ctorsion_t& c = ctorsions[torsions[i].code - 1];
        ids.push_back({torsions[i].ai - 1, torsions[i].aj - 1,
                       torsions[i].ak - 1, torsions[i].al - 1});
        params.push_back({static_cast<real_t>(c.k),
                          static_cast<real_t>(c.n),
                          static_cast<real_t>(deg2rad(c.d)),
                          static_cast<real_t>(c.paths)});  // {k, n, d(rad), paths}
        eslot.push_back(i < ctx.n_torsions_solute ? bonded_p_slot() : bonded_w_slot());
    }

    data_.torsion.n = static_cast<int>(ids.size());
    data_.torsion.ids = HostDeviceBuffer<dihe_idx_t>::from_vector(ids, run_gpu);
    data_.torsion.params = HostDeviceBuffer<torsion_param_t>::from_vector(params, run_gpu);
    data_.torsion.eslot = HostDeviceBuffer<int>::from_vector(eslot, run_gpu);
}
void BondedForce::build_impropers(Context& ctx) {
    const bool run_gpu = ctx.command_info.requested_gpu;
    auto* impropers = ctx.impropers->cpu_data_p;
    auto* cimpropers = ctx.cimpropers->cpu_data_p;

    std::vector<dihe_idx_t> ids;
    std::vector<dparam2_t> params;
    std::vector<int> eslot;
    ids.reserve(ctx.n_impropers);
    params.reserve(ctx.n_impropers);
    eslot.reserve(ctx.n_impropers);

    for (int i = 0; i < ctx.n_impropers; i++) {
        const cimproper_t& c = cimpropers[impropers[i].code - 1];
        ids.push_back({impropers[i].ai - 1, impropers[i].aj - 1,
                       impropers[i].ak - 1, impropers[i].al - 1});
        params.push_back({c.k, deg2rad(c.phi0)});  // {k, phi0(rad)}
        eslot.push_back(i < ctx.n_impropers_solute ? bonded_p_slot() : bonded_w_slot());
    }

    data_.improper.n = static_cast<int>(ids.size());
    data_.improper.ids = HostDeviceBuffer<dihe_idx_t>::from_vector(ids, run_gpu);
    data_.improper.params = HostDeviceBuffer<dparam2_t>::from_vector(params, run_gpu);
    data_.improper.eslot = HostDeviceBuffer<int>::from_vector(eslot, run_gpu);
}
