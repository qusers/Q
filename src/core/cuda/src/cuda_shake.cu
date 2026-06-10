#include <map>
#include <unordered_map>

#include "constants.h"
#include "cuda_runtime_utility.h"
#include "cuda_shake.cuh"

namespace {
const int kShakeThreads = 64;

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
    const real_t height = std::sqrt(roh2 - 0.25 * rhh2);
    const real_t total_mass = wo + 2.0 * wh;
    const real_t com_y = 2.0 * wh * height / total_mass;

    return ShakeFastWater{
        o, h1, h2, com_y, 1.0 / com_y, height - com_y, 0.5 * rhh, rhh, rhh * rhh, wo / total_mass, wh / total_mass};
}

__device__ real_t clamp_nonnegative(real_t value) {
    return value > static_cast<real_t>(0) ? value : static_cast<real_t>(0);
}

}  // namespace

void CudaShake::find_shake_fast_water(Context& ctx, std::vector<bool>& optimized) {
    std::vector<ShakeFastWater> shake_fast_waters;
    const auto* winv = ctx.winv->cpu_data_p;

    std::map<BondKey, int> shake_bonds_lookup_table;

    const auto& shake_bond = data_.shake_bonds->cpu_data_p;
    for (int i = 0; i < data_.n_constraints; i++) {
        int ai = shake_bond[i].ai - 1;
        int aj = shake_bond[i].aj - 1;

        shake_bonds_lookup_table[make_bond_key(ai, aj)] = i;
    }

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
    this->shake_fast_waters = HostDeviceBuffer<ShakeFastWater>::from_vector(shake_fast_waters, true);
}

void CudaShake::find_shake_network(Context& ctx, std::vector<bool>& optimized) {
    std::vector<std::vector<int>> atom_to_bonds(ctx.n_atoms);

    auto* shake_bonds = data_.shake_bonds->cpu_data_p;
    auto* heavy = ctx.heavy->cpu_data_p;
    auto* winv = ctx.winv->cpu_data_p;

    for (int i = 0; i < data_.n_constraints; i++) {
        if (optimized[i]) continue;
        atom_to_bonds[shake_bonds[i].ai - 1].push_back(i);
        atom_to_bonds[shake_bonds[i].aj - 1].push_back(i);
    }

    std::vector<ShakeNetwork> shake_networks;
    std::vector<bool> bond_visited(data_.n_constraints, false);
    std::vector<bool> atom_visited(ctx.n_atoms, false);

    for (int i = 0; i < data_.n_constraints; i++) {
        if (optimized[i] || bond_visited[i]) continue;

        std::vector<int> component_bonds, component_atoms;
        std::vector<int> stk = {i};
        bond_visited[i] = true;

        while (!stk.empty()) {
            const int bi = stk.back();
            stk.pop_back();

            component_bonds.push_back(bi);

            const int atoms[2] = {shake_bonds[bi].ai - 1, shake_bonds[bi].aj - 1};

            for (int atom : atoms) {
                if (atom_visited[atom]) continue;
                atom_visited[atom] = true;
                component_atoms.push_back(atom);

                for (int next_bond : atom_to_bonds[atom]) {
                    if (!bond_visited[next_bond]) {
                        bond_visited[next_bond] = true;
                        stk.push_back(next_bond);
                    }
                }
            }
        }

        if (component_bonds.size() > 3 || component_atoms.size() != component_bonds.size() + 1) continue;

        int center = -1;

        for (int atom : component_atoms) {
            if (heavy[atom]) {
                if (center != -1) {
                    center = -2;
                    break;
                }
                center = atom;
            }
        }
        if (center < 0) continue;

        ShakeNetwork network = {};
        network.center = center;
        network.center_winv = winv[center];
        bool valid = true;
        for (int bi : component_bonds) {
            const int ai = shake_bonds[bi].ai - 1;
            const int aj = shake_bonds[bi].aj - 1;
            int hydrogen = -1;
            if (ai == center && !heavy[aj]) {
                hydrogen = aj;
            } else if (aj == center && !heavy[ai]) {
                hydrogen = ai;
            } else {
                valid = false;
                break;
            }

            int hidx = network.n_hydrogens;
            if (hidx >= 3) {
                valid = false;
                break;
            }

            network.hydrogens[hidx] = hydrogen;
            network.dist2[hidx] = shake_bonds[bi].dist2;
            network.hydrogen_winv[hidx] = winv[hydrogen];
            network.n_hydrogens++;
        }

        if (!valid || network.n_hydrogens == 0) continue;
        shake_networks.push_back(network);
        for (int bi : component_bonds) optimized[bi] = true;
    }
    this->shake_networks = HostDeviceBuffer<ShakeNetwork>::from_vector(shake_networks, true);
}

void CudaShake::find_fallback_shake_bond(Context& ctx, std::vector<bool>& optimized) {
    std::vector<std::vector<int>> atom_to_bonds(ctx.n_atoms);
    auto* shake_bonds = data_.shake_bonds->cpu_data_p;
    for (int i = 0; i < data_.n_constraints; i++) {
        if (optimized[i]) continue;
        atom_to_bonds[shake_bonds[i].ai - 1].push_back(i);
        atom_to_bonds[shake_bonds[i].aj - 1].push_back(i);
    }

    std::vector<std::vector<ShakeBond>> fallback_bonds_by_color;
    std::vector<bool> visited(data_.n_constraints, false);

    for (int i = 0; i < data_.n_constraints; i++) {
        if (optimized[i] || visited[i]) continue;

        std::vector<int> component_bonds;
        std::vector<int> stk = {i};
        visited[i] = true;

        while (!stk.empty()) {
            const int bi = stk.back();
            stk.pop_back();
            component_bonds.push_back(bi);

            const int atoms[2] = {shake_bonds[bi].ai - 1, shake_bonds[bi].aj - 1};

            for (int atom : atoms) {
                for (int next_bond : atom_to_bonds[atom]) {
                    if (!visited[next_bond]) {
                        visited[next_bond] = true;
                        stk.push_back(next_bond);
                    }
                }
            }
        }

        std::unordered_map<int, std::vector<int>> atom_to_colors;
        std::vector<int> component_colors;
        int component_n_colors = 0;
        for (int bi : component_bonds) {
            std::vector<bool> used(component_n_colors + 1, false);

            const int ai = shake_bonds[bi].ai - 1;
            const int aj = shake_bonds[bi].aj - 1;

            auto mark_used = [&](int atom) {
                auto it = atom_to_colors.find(atom);
                if (it == atom_to_colors.end()) return;
                for (int color : it->second) used[color] = 1;
            };

            mark_used(ai);
            mark_used(aj);

            int color = 0;
            while (color < component_n_colors && used[color]) color++;
            if (color == component_n_colors) component_n_colors++;

            component_colors.push_back(color);
            atom_to_colors[ai].push_back(color);
            atom_to_colors[aj].push_back(color);
        }

        if (component_n_colors > fallback_bonds_by_color.size()) {
            fallback_bonds_by_color.resize(component_n_colors);
        }

        for (int i = 0; i < component_bonds.size(); i++) {
            const int bi = component_bonds[i];
            fallback_bonds_by_color[component_colors[i]].push_back(shake_bonds[bi]);
        }
    }

    std::vector<int> fallback_color_offsets(fallback_bonds_by_color.size() + 1, 0);
    std::vector<ShakeBond> fallback_shake_bonds;
    for (int color = 0; color < fallback_bonds_by_color.size(); color++) {
        fallback_color_offsets[color] = fallback_shake_bonds.size();
        fallback_shake_bonds.insert(
            fallback_shake_bonds.end(),
            fallback_bonds_by_color[color].begin(),
            fallback_bonds_by_color[color].end());
    }
    fallback_color_offsets[fallback_bonds_by_color.size()] = fallback_shake_bonds.size();

    this->fallback_color_offsets = fallback_color_offsets;
    this->fallback_shake_bonds = HostDeviceBuffer<ShakeBond>::from_vector(fallback_shake_bonds, true);
}

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

    std::vector<bool> optimized(data_.n_constraints, false);

    find_shake_fast_water(ctx, optimized);
    find_shake_network(ctx, optimized);
    find_fallback_shake_bond(ctx, optimized);
    fallback_unconverged = std::make_unique<HostDeviceBuffer<int>>(1, true, true);
    is_init_backend = true;
}

__global__ void calc_fast_water_shake_kernel(
    int n_fast_waters,
    ShakeFastWater* fast_waters,
    coord_t* coords,
    coord_t* xcoords) {
    /*
    todo: Need to think...
    */
    int idx = blockIdx.x * blockDim.x + threadIdx.x;
    if (idx >= n_fast_waters) return;

    auto water = fast_waters[idx];
    const int o = water.o;
    const int h1 = water.h1;
    const int h2 = water.h2;

    const real_t xo0 = xcoords[o].x;
    const real_t yo0 = xcoords[o].y;
    const real_t zo0 = xcoords[o].z;
    const real_t xh10 = xcoords[h1].x;
    const real_t yh10 = xcoords[h1].y;
    const real_t zh10 = xcoords[h1].z;
    const real_t xh20 = xcoords[h2].x;
    const real_t yh20 = xcoords[h2].y;
    const real_t zh20 = xcoords[h2].z;

    const real_t xo1 = coords[o].x;
    const real_t yo1 = coords[o].y;
    const real_t zo1 = coords[o].z;
    const real_t xh11 = coords[h1].x;
    const real_t yh11 = coords[h1].y;
    const real_t zh11 = coords[h1].z;
    const real_t xh21 = coords[h2].x;
    const real_t yh21 = coords[h2].y;
    const real_t zh21 = coords[h2].z;

    const real_t xb0 = xh10 - xo0;
    const real_t yb0 = yh10 - yo0;
    const real_t zb0 = zh10 - zo0;
    const real_t xc0 = xh20 - xo0;
    const real_t yc0 = yh20 - yo0;
    const real_t zc0 = zh20 - zo0;

    const real_t xcom = xo1 * water.wo_div_wohh + (xh11 + xh21) * water.wh_div_wohh;
    const real_t ycom = yo1 * water.wo_div_wohh + (yh11 + yh21) * water.wh_div_wohh;
    const real_t zcom = zo1 * water.wo_div_wohh + (zh11 + zh21) * water.wh_div_wohh;

    const real_t xa1 = xo1 - xcom;
    const real_t ya1 = yo1 - ycom;
    const real_t za1 = zo1 - zcom;
    const real_t xb1 = xh11 - xcom;
    const real_t yb1 = yh11 - ycom;
    const real_t zb1 = zh11 - zcom;
    const real_t xc1 = xh21 - xcom;
    const real_t yc1 = yh21 - ycom;
    const real_t zc1 = zh21 - zcom;

    const real_t xakszd = yb0 * zc0 - zb0 * yc0;
    const real_t yakszd = zb0 * xc0 - xb0 * zc0;
    const real_t zakszd = xb0 * yc0 - yb0 * xc0;
    const real_t xaksxd = ya1 * zakszd - za1 * yakszd;
    const real_t yaksxd = za1 * xakszd - xa1 * zakszd;
    const real_t zaksxd = xa1 * yakszd - ya1 * xakszd;
    const real_t xaksyd = yakszd * zaksxd - zakszd * yaksxd;
    const real_t yaksyd = zakszd * xaksxd - xakszd * zaksxd;
    const real_t zaksyd = xakszd * yaksxd - yakszd * xaksxd;

    const real_t ax2 = xaksxd * xaksxd + yaksxd * yaksxd + zaksxd * zaksxd;
    const real_t ay2 = xaksyd * xaksyd + yaksyd * yaksyd + zaksyd * zaksyd;
    const real_t az2 = xakszd * xakszd + yakszd * yakszd + zakszd * zakszd;
    if (ax2 <= static_cast<real_t>(0) || ay2 <= static_cast<real_t>(0) || az2 <= static_cast<real_t>(0)) return;

    const real_t axlng_inv = rsqrt(ax2);
    const real_t aylng_inv = rsqrt(ay2);
    const real_t azlng_inv = rsqrt(az2);

    const real_t trns11 = xaksxd * axlng_inv;
    const real_t trns21 = yaksxd * axlng_inv;
    const real_t trns31 = zaksxd * axlng_inv;
    const real_t trns12 = xaksyd * aylng_inv;
    const real_t trns22 = yaksyd * aylng_inv;
    const real_t trns32 = zaksyd * aylng_inv;
    const real_t trns13 = xakszd * azlng_inv;
    const real_t trns23 = yakszd * azlng_inv;
    const real_t trns33 = zakszd * azlng_inv;

    const real_t xb0d = trns11 * xb0 + trns21 * yb0 + trns31 * zb0;
    const real_t yb0d = trns12 * xb0 + trns22 * yb0 + trns32 * zb0;
    const real_t xc0d = trns11 * xc0 + trns21 * yc0 + trns31 * zc0;
    const real_t yc0d = trns12 * xc0 + trns22 * yc0 + trns32 * zc0;
    const real_t za1d = trns13 * xa1 + trns23 * ya1 + trns33 * za1;
    const real_t xb1d = trns11 * xb1 + trns21 * yb1 + trns31 * zb1;
    const real_t yb1d = trns12 * xb1 + trns22 * yb1 + trns32 * zb1;
    const real_t zb1d = trns13 * xb1 + trns23 * yb1 + trns33 * zb1;
    const real_t xc1d = trns11 * xc1 + trns21 * yc1 + trns31 * zc1;
    const real_t yc1d = trns12 * xc1 + trns22 * yc1 + trns32 * zc1;
    const real_t zc1d = trns13 * xc1 + trns23 * yc1 + trns33 * zc1;

    const real_t sinphi = za1d * water.ra_inv;
    const real_t cosphi = sqrt(clamp_nonnegative(static_cast<real_t>(1) - sinphi * sinphi));
    if (cosphi <= static_cast<real_t>(0)) return;
    const real_t sinpsi = (zb1d - zc1d) / (water.rhh * cosphi);
    const real_t cospsi = sqrt(clamp_nonnegative(static_cast<real_t>(1) - sinpsi * sinpsi));

    const real_t ya2d = water.ra * cosphi;
    real_t xb2d = -water.rc * cospsi;
    const real_t yb2d = -water.rb * cosphi - water.rc * sinpsi * sinphi;
    const real_t yc2d = -water.rb * cosphi + water.rc * sinpsi * sinphi;
    xb2d = -static_cast<real_t>(0.5) *
           sqrt(clamp_nonnegative(water.rhh2 - (yb2d - yc2d) * (yb2d - yc2d) - (zb1d - zc1d) * (zb1d - zc1d)));

    const real_t alpa = xb2d * (xb0d - xc0d) + yb0d * yb2d + yc0d * yc2d;
    const real_t beta = xb2d * (yc0d - yb0d) + xb0d * yb2d + xc0d * yc2d;
    const real_t gama = xb0d * yb1d - xb1d * yb0d + xc0d * yc1d - xc1d * yc0d;
    const real_t al2be2 = alpa * alpa + beta * beta;
    if (al2be2 <= static_cast<real_t>(0)) return;

    const real_t sinthe = (alpa * gama - beta * sqrt(clamp_nonnegative(al2be2 - gama * gama))) / al2be2;
    const real_t costhe = sqrt(clamp_nonnegative(static_cast<real_t>(1) - sinthe * sinthe));

    const real_t xa3d = -ya2d * sinthe;
    const real_t ya3d = ya2d * costhe;
    const real_t za3d = za1d;
    const real_t xb3d = xb2d * costhe - yb2d * sinthe;
    const real_t yb3d = xb2d * sinthe + yb2d * costhe;
    const real_t zb3d = zb1d;
    const real_t xc3d = -xb2d * costhe - yc2d * sinthe;
    const real_t yc3d = -xb2d * sinthe + yc2d * costhe;
    const real_t zc3d = zc1d;

    coords[o].x = xcom + trns11 * xa3d + trns12 * ya3d + trns13 * za3d;
    coords[o].y = ycom + trns21 * xa3d + trns22 * ya3d + trns23 * za3d;
    coords[o].z = zcom + trns31 * xa3d + trns32 * ya3d + trns33 * za3d;
    coords[h1].x = xcom + trns11 * xb3d + trns12 * yb3d + trns13 * zb3d;
    coords[h1].y = ycom + trns21 * xb3d + trns22 * yb3d + trns23 * zb3d;
    coords[h1].z = zcom + trns31 * xb3d + trns32 * yb3d + trns33 * zb3d;
    coords[h2].x = xcom + trns11 * xc3d + trns12 * yc3d + trns13 * zc3d;
    coords[h2].y = ycom + trns21 * xc3d + trns22 * yc3d + trns23 * zc3d;
    coords[h2].z = zcom + trns31 * xc3d + trns32 * yc3d + trns33 * zc3d;

}

__global__ void calc_h_star_shake_kernel(
    int n_shake_networks,
    ShakeNetwork* networks,
    coord_t* coords,
    coord_t* xcoords) {
    int idx = blockIdx.x * blockDim.x + threadIdx.x;
    if (idx >= n_shake_networks) return;

    const auto network = networks[idx];
    coord_t center = coords[network.center];
    coord_t hydrogens[3];
    coord_t old_vectors[3];

    for (int i = 0; i < network.n_hydrogens; i++) {
        const int h = network.hydrogens[i];
        hydrogens[i] = coords[h];
        old_vectors[i].x = xcoords[network.center].x - xcoords[h].x;
        old_vectors[i].y = xcoords[network.center].y - xcoords[h].y;
        old_vectors[i].z = xcoords[network.center].z - xcoords[h].z;
    }

    bool converged = false;
    int n_iterations = 0;
    do {
        converged = true;
        for (int i = 0; i < network.n_hydrogens; i++) {
            coord_t current_vector;
            current_vector.x = center.x - hydrogens[i].x;
            current_vector.y = center.y - hydrogens[i].y;
            current_vector.z = center.z - hydrogens[i].z;
            const real_t current_dist2 = current_vector.x * current_vector.x +
                                         current_vector.y * current_vector.y +
                                         current_vector.z * current_vector.z;
            const real_t diff = network.dist2[i] - current_dist2;
            if (fabs(diff) < shake_tol * network.dist2[i]) continue;

            converged = false;
            const real_t scp = current_vector.x * old_vectors[i].x +
                               current_vector.y * old_vectors[i].y +
                               current_vector.z * old_vectors[i].z;
            const real_t inv_mass_sum = network.center_winv + network.hydrogen_winv[i];
            if (scp <= network.dist2[i] * static_cast<real_t>(1.0e-6) || inv_mass_sum == static_cast<real_t>(0)) continue;

            const real_t corr = diff / (static_cast<real_t>(2) * scp * inv_mass_sum);
            const real_t center_scale = corr * network.center_winv;
            const real_t hydrogen_scale = corr * network.hydrogen_winv[i];
            center.x += old_vectors[i].x * center_scale;
            center.y += old_vectors[i].y * center_scale;
            center.z += old_vectors[i].z * center_scale;
            hydrogens[i].x -= old_vectors[i].x * hydrogen_scale;
            hydrogens[i].y -= old_vectors[i].y * hydrogen_scale;
            hydrogens[i].z -= old_vectors[i].z * hydrogen_scale;
        }
        n_iterations++;
    } while (n_iterations < shake_max_iter && !converged);

    if (!converged) {
        for (int i = 0; i < network.n_hydrogens; i++) {
            const real_t dx = center.x - hydrogens[i].x;
            const real_t dy = center.y - hydrogens[i].y;
            const real_t dz = center.z - hydrogens[i].z;
            const real_t dist2 = dx * dx + dy * dy + dz * dz;
            printf(">>> Shake failed, i = %d,j = %d, d = %f, d0 = %f\n",
                   network.center,
                   network.hydrogens[i],
                   static_cast<double>(sqrt(dist2)),
                   static_cast<double>(sqrt(network.dist2[i])));
        }
        return;
    }

    coords[network.center] = center;
    for (int i = 0; i < network.n_hydrogens; i++) {
        coords[network.hydrogens[i]] = hydrogens[i];
    }
}

__global__ void print_fallback_shake_failures_kernel(
    int n_shakes,
    ShakeBond* shake_bonds,
    coord_t* coords) {
    const int idx = blockIdx.x * blockDim.x + threadIdx.x;
    if (idx >= n_shakes) return;

    const int ai = shake_bonds[idx].ai - 1;
    const int aj = shake_bonds[idx].aj - 1;
    const real_t dx = coords[ai].x - coords[aj].x;
    const real_t dy = coords[ai].y - coords[aj].y;
    const real_t dz = coords[ai].z - coords[aj].z;
    const real_t dist2 = dx * dx + dy * dy + dz * dz;
    if (fabs(shake_bonds[idx].dist2 - dist2) >= shake_tol * shake_bonds[idx].dist2) {
        printf(">>> Shake failed, i = %d,j = %d, d = %f, d0 = %f\n",
               ai,
               aj,
               static_cast<double>(sqrt(dist2)),
               static_cast<double>(sqrt(shake_bonds[idx].dist2)));
    }
}

__global__ void calc_fallback_shake_color_kernel(
    int n_shakes,
    ShakeBond* shake_bonds,
    coord_t* coords,
    coord_t* xcoords,
    real_t* winv,
    int* unconverged) {
    const int idx = blockIdx.x * blockDim.x + threadIdx.x;
    if (idx >= n_shakes) return;

    ShakeBond& shake_bond = shake_bonds[idx];
    const int ai = shake_bond.ai - 1;
    const int aj = shake_bond.aj - 1;
    const real_t xij_x = coords[ai].x - coords[aj].x;
    const real_t xij_y = coords[ai].y - coords[aj].y;
    const real_t xij_z = coords[ai].z - coords[aj].z;
    const real_t xij2 = xij_x * xij_x + xij_y * xij_y + xij_z * xij_z;
    const real_t diff = shake_bond.dist2 - xij2;
    if (fabs(diff) < shake_tol * shake_bond.dist2) return;

    atomicExch(unconverged, 1);
    const real_t xxij_x = xcoords[ai].x - xcoords[aj].x;
    const real_t xxij_y = xcoords[ai].y - xcoords[aj].y;
    const real_t xxij_z = xcoords[ai].z - xcoords[aj].z;
    const real_t scp = xij_x * xxij_x + xij_y * xxij_y + xij_z * xxij_z;
    const real_t inv_mass_sum = winv[ai] + winv[aj];
    if (scp == static_cast<real_t>(0) || inv_mass_sum == static_cast<real_t>(0)) return;

    const real_t corr = diff / (static_cast<real_t>(2) * scp * inv_mass_sum);
    const real_t ai_scale = corr * winv[ai];
    const real_t aj_scale = corr * winv[aj];
    coords[ai].x += xxij_x * ai_scale;
    coords[ai].y += xxij_y * ai_scale;
    coords[ai].z += xxij_z * ai_scale;
    coords[aj].x -= xxij_x * aj_scale;
    coords[aj].y -= xxij_y * aj_scale;
    coords[aj].z -= xxij_z * aj_scale;
}

void CudaShake::apply_to(Context& ctx, coord_t* d_coords, coord_t* d_xcoords) {
    if (shake_fast_waters->length > 0) {
        const int grid_blocks = (shake_fast_waters->length + kShakeThreads - 1) / kShakeThreads;
        calc_fast_water_shake_kernel<<<grid_blocks, kShakeThreads>>>(shake_fast_waters->length, shake_fast_waters->gpu_data_p, d_coords, d_xcoords);
    }
    if (shake_networks->length > 0) {
        const int grid_blocks = (shake_networks->length + kShakeThreads - 1) / kShakeThreads;
        calc_h_star_shake_kernel<<<grid_blocks, kShakeThreads>>>(shake_networks->length, shake_networks->gpu_data_p, d_coords, d_xcoords);
    }
    if (fallback_shake_bonds->length > 0) {
        int color_groups_num = fallback_color_offsets.size() - 1;

        for (int n_iterations = 0; n_iterations < shake_max_iter; n_iterations++) {
            fallback_unconverged->cpu_data_p[0] = false;
            fallback_unconverged->upload();
            for (int i = 0; i < color_groups_num; i++) {
                const int offset = fallback_color_offsets[i];
                const int bond_num = fallback_color_offsets[i + 1] - offset;
                const int grid_blocks = (bond_num + kShakeThreads - 1) / kShakeThreads;
                calc_fallback_shake_color_kernel<<<grid_blocks, kShakeThreads>>>(
                    bond_num,
                    fallback_shake_bonds->gpu_data_p + offset,
                    d_coords,
                    d_xcoords,
                    ctx.winv->gpu_data_p,
                    fallback_unconverged->gpu_data_p);
            }
            fallback_unconverged->download();
            if (fallback_unconverged->cpu_data_p[0] == false) {
                break;
            }
        }

        if (fallback_unconverged->cpu_data_p[0]) {
            int n_fallback_constraints = fallback_shake_bonds->length;
            const int grid_blocks = (n_fallback_constraints + kShakeThreads - 1) / kShakeThreads;
            print_fallback_shake_failures_kernel<<<grid_blocks, kShakeThreads>>>(
                n_fallback_constraints,
                fallback_shake_bonds->gpu_data_p,
                d_coords);
        }
    }
}

void CudaShake::cleanup() {
}
