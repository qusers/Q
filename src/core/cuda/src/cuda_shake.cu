#include <map>

#include "cuda_runtime_utility.h"
#include "cuda_shake.cuh"

namespace {
using BondKey = std::pair<int, int>;

BondKey make_bond_key(int ai, int aj) {
    if (ai > aj) std::swap(ai, aj);
    return {ai, aj};
}

int find_bond_index(const std::map<BondKey, int>& bond_by_atoms, int ai, int aj) {
    auto it = bond_by_atoms.find(make_bond_key(ai, aj));
    return it == bond_by_atoms.end() ? -1 : it->second;
}

ShakeFastWater make_fast_water(int o, int h1, int h2, real_t roh2, real_t rhh2, real_t wo, real_t wh) {
    const real_t roh = std::sqrt(roh2);
    const real_t rhh = std::sqrt(rhh2);
    const real_t height = std::sqrt(roh2 - 0.25 * rhh);
    const real_t total_mass = wo + 2.0 * wh;
    const real_t com_y = -2.0 * wh * height / total_mass;

    return ShakeFastWater{
        o, h1, h2, 
        
    }


}

}  // namespace

void CudaShake::apply(Context& ctx) {
    auto* coords = ctx.coords->gpu_data_p;
    auto* xcoords = ctx.xcoords->gpu_data_p;
    apply_to(ctx, coords, xcoords);
}

void CudaShake::initial_shake(Context& ctx) {
    auto* d_coords = ctx.coords->gpu_data_p;
    auto* d_xcoords = ctx.xcoords->gpu_data_p;

    check_cuda(cudaMemcpy(d_xcoords, d_coords, ctx.n_atoms * sizeof(coord_t), cudaMemcpyDeviceToDevice));
    apply_to(ctx, d_coords, d_xcoords);

    ctx.coords->download();

    auto* coords = ctx.coords->cpu_data_p;
    auto* velocities = ctx.velocities->cpu_data_p;
    auto* xcoords = ctx.xcoords->cpu_data_p;
    for (int i = 0; i < ctx.n_atoms; i++) {
        xcoords[i].x = coords[i].x - ctx.dt * velocities[i].x;
        xcoords[i].y = coords[i].y - ctx.dt * velocities[i].y;
        xcoords[i].z = coords[i].z - ctx.dt * velocities[i].z;
    }

    ctx.xcoords->upload();

    apply_to(ctx, d_xcoords, d_coords);

    ctx.xcoords->download();

    for (int i = 0; i < ctx.n_atoms; i++) {
        velocities[i].x = (coords[i].x - xcoords[i].x) / ctx.dt;
        velocities[i].y = (coords[i].y - xcoords[i].y) / ctx.dt;
        velocities[i].z = (coords[i].z - xcoords[i].z) / ctx.dt;
    }
    ctx.velocities->upload();
}

void CudaShake::init_backend(Context& ctx) {
    if (is_init_backend) {
        return;
    }

    std::map<BondKey, int> shake_bonds_lookup_table;
    std::vector<int> atom_degree(ctx.n_atoms);

    const auto& shake_bond = data_.shake_bonds->cpu_data_p;
    for (int i = 0; i < data_.n_constraints; i++) {
        int ai = shake_bond[i].ai - 1;
        int aj = shake_bond[i].aj - 1;

        shake_bonds_lookup_table[make_bond_key(ai, aj)] = i;
        atom_degree[ai]++;
        atom_degree[aj]++;
    }

    /*
    for fast water
    */

    std::vector<ShakeFastWater> shake_fast_waters;
    const auto* winv = ctx.winv->cpu_data_p;
    std::vector<bool> optimized(data_.n_constraints, false);

    for (int w = 0; w < ctx.n_waters; w++) {
        const int o = ctx.n_atoms_solute + 3 * w;
        const int h1 = o + 1;
        const int h2 = o + 2;

        const int oh1 = find_bond_index(shake_bonds_lookup_table, o, h1);
        const int oh2 = find_bond_index(shake_bonds_lookup_table, o, h2);
        const int hh = find_bond_index(shake_bonds_lookup_table, h1, h2);

        if (oh1 < 0 || oh2 < 0 || hh < 0) continue;
        if (shake_bond[oh1].dist2 != shake_bond[oh2].dist2) continue;
        // triangle height from oxygen to the H-H midpoint height = sqrt(rOH * rOH - 0.25 * rHH * rHH);
        if (shake_bond[oh1].dist2 <= 0.25 * shake_bond[hh].dist2) continue;

        const real_t wo = 1.0 / winv[o];
        const real_t wh1 = 1.0 / winv[h1];
        const real_t wh2 = 1.0 / winv[h2];
        if (wh1 != wh2) continue;
        shake_fast_waters.push_back(make_fast_water(o, h1, h2, shake_bond[oh1].dist2, shake_bond[hh].dist2, wo, wh1));
        optimized[oh1] = true;
        optimized[oh2] = true;
        optimized[hh] = true;
    }
}

void CudaShake::apply_to(Context& ctx, coord_t* d_coords, coord_t* d_xcoords) {
}

void CudaShake::cleanup() {
}
