#include "restraint_force.h"

#include "constants.h"

void RestraintForce::init(Context& ctx, const std::vector<restrpos_t>& restrspos,
                          const std::vector<restrseq_t>& restrseqs,
                          const std::vector<restrdis_t>& restrdists,
                          const std::vector<restrang_t>& restrangs,
                          const std::vector<restrwall_t>& restrwalls) {
    build_posrestr(ctx);
    build_restrpos(ctx, restrspos);
    build_restrseq(ctx, restrseqs);
    build_restrdis(ctx, restrdists);
    build_restrang(ctx, restrangs);
    build_restrwall(ctx, restrwalls);
    init_backend(ctx);
}

// pshell + fixed atoms: derived (SoA), harmonic to coords_init, isotropic k, no FEP.
void RestraintForce::build_posrestr(Context& ctx) {
    const bool run_gpu = ctx.command_info.requested_gpu;
    const auto* shell = ctx.shell->cpu_data_p;
    const auto* excluded = ctx.excluded->cpu_data_p;

    std::vector<int> ids, eslot_a, eslot_b;
    std::vector<double> k;
    for (int i = 0; i < ctx.n_atoms_solute; i++) {
        if (!(shell[i] || excluded[i])) continue;  // atom should outside boundary shell or frozen
        ids.push_back(i);
        k.push_back(excluded[i] ? k_fix : k_pshell);
        eslot_a.push_back(excluded[i] ? E_RESTR_FIX : -1);
        eslot_b.push_back(shell[i] ? E_RESTR_SHELL : -1);
    }
    auto& p = data_.posrestr;
    p.n = static_cast<int>(ids.size());
    p.ids = HostDeviceBuffer<int>::from_vector(ids, run_gpu);
    p.k = HostDeviceBuffer<double>::from_vector(k, run_gpu);
    p.eslot_a = HostDeviceBuffer<int>::from_vector(eslot_a, run_gpu);
    p.eslot_b = HostDeviceBuffer<int>::from_vector(eslot_b, run_gpu);
}

void RestraintForce::build_restrpos(Context& ctx, const std::vector<restrpos_t>& restrspos) {
    const bool run_gpu = ctx.command_info.requested_gpu;
    data_.restrpos.n = restrspos.size();
    if (data_.restrpos.n == 0) return;
    data_.restrpos.recs = HostDeviceBuffer<restrpos_t>::from_vector(restrspos, run_gpu);
}

void RestraintForce::build_restrseq(Context& ctx, const std::vector<restrseq_t>& restrseqs) {
    const bool run_gpu = ctx.command_info.requested_gpu;
    data_.restrseq.n = restrseqs.size();
    if (data_.restrseq.n == 0) return;
    data_.restrseq.recs = HostDeviceBuffer<restrseq_t>::from_vector(restrseqs, run_gpu);
}

void RestraintForce::build_restrdis(Context& ctx, const std::vector<restrdis_t>& restrdists) {
    const bool run_gpu = ctx.command_info.requested_gpu;
    data_.restrdis.n = restrdists.size();
    if (data_.restrdis.n == 0) return;
    data_.restrdis.recs = HostDeviceBuffer<restrdis_t>::from_vector(restrdists, run_gpu);
}

void RestraintForce::build_restrang(Context& ctx, const std::vector<restrang_t>& restrangs) {
    const bool run_gpu = ctx.command_info.requested_gpu;
    data_.restrang.n = restrangs.size();
    if (data_.restrang.n == 0) return;
    data_.restrang.recs = HostDeviceBuffer<restrang_t>::from_vector(restrangs, run_gpu);
}

void RestraintForce::build_restrwall(Context& ctx, const std::vector<restrwall_t>& restrwalls) {
    const bool run_gpu = ctx.command_info.requested_gpu;
    data_.restrwall.n = restrwalls.size();
    if (data_.restrwall.n == 0) return;
    data_.restrwall.recs = HostDeviceBuffer<restrwall_t>::from_vector(restrwalls, run_gpu);
}
