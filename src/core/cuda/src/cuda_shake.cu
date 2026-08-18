#include <cooperative_groups.h>

#include <algorithm>
#include <cstdio>
#include <cstdlib>
#include <map>
#include <unordered_map>

#include "constants.h"
#include "cuda_runtime_utility.h"
#include "cuda_shake.cuh"
namespace cg = cooperative_groups;

namespace {

constexpr double kMinShakeProjectionFactor = 1.0e-6;
constexpr double kMaxShakeCorrectionScale = 1;
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

ShakeFastWater make_fast_water(int o, int h1, int h2, double roh2, double rhh2, double wo, double wh) {
    const double roh = sqrt(roh2);
    const double rhh = sqrt(rhh2);
    const double height = sqrt(roh2 - 0.25 * rhh2);
    const double total_mass = wo + 2.0 * wh;
    const double com_y = 2.0 * wh * height / total_mass;

    return ShakeFastWater{
        o, h1, h2, com_y, 1.0 / com_y, height - com_y, 0.5 * rhh, rhh, rhh * rhh, wo / total_mass, wh / total_mass};
}

__device__ double clamp_nonnegative(double value) {
    return value > 0.0 ? value : 0.0;
}

}  // namespace

void CudaShake::find_serial_q_molecule_bonds(Context& ctx, std::vector<bool>& optimized) {
    std::vector<bool> serial_molecules(ctx.n_molecules(), false);
    for (const int atom : ctx.q_atoms) {
        const int atom_number = atom + 1;
        auto next = std::upper_bound(ctx.molecules.begin(), ctx.molecules.end(), atom_number);
        if (next != ctx.molecules.begin()) {
            serial_molecules[static_cast<int>(next - ctx.molecules.begin()) - 1] = true;
        }
    }

    std::vector<int> molecule_constraint_counts(ctx.n_molecules(), 0);
    std::vector<ConstraintBond> serial_bonds;
    const ConstraintBond* shake_bonds = data_.constraint_bonds->cpu_data_p;
    int current_mol = 0;
    for (int i = 0; i < data_.n_constraints; i++) {
        const int atom_number = shake_bonds[i].ai;
        while (current_mol + 1 < ctx.n_molecules() &&
               atom_number >= ctx.molecules[current_mol + 1]) {
            current_mol++;
        }
        if (!serial_molecules[current_mol]) continue;

        optimized[i] = true;
        molecule_constraint_counts[current_mol]++;
        serial_bonds.push_back(shake_bonds[i]);
    }

    serial_q_solver_.init(molecule_constraint_counts, serial_bonds);
}

void CudaShake::find_shake_fast_water(Context& ctx, std::vector<bool>& optimized) {
    std::vector<ShakeFastWater> shake_fast_waters;
    const auto* winv = ctx.winv->cpu_data_p;

    std::map<BondKey, int> shake_bonds_lookup_table;

    const auto& shake_bond = data_.constraint_bonds->cpu_data_p;
    for (int i = 0; i < data_.n_constraints; i++) {
        int ai = shake_bond[i].ai - 1;
        int aj = shake_bond[i].aj - 1;

        shake_bonds_lookup_table[make_bond_key(ai, aj)] = i;
    }

    for (int w = 0; w < ctx.n_waters(); w++) {
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

        const double wo = 1.0 / winv[o];
        const double wh1 = 1.0 / winv[h1];
        const double wh2 = 1.0 / winv[h2];
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

    auto* shake_bonds = data_.constraint_bonds->cpu_data_p;
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
    auto* shake_bonds = data_.constraint_bonds->cpu_data_p;
    for (int i = 0; i < data_.n_constraints; i++) {
        if (optimized[i]) continue;
        atom_to_bonds[shake_bonds[i].ai - 1].push_back(i);
        atom_to_bonds[shake_bonds[i].aj - 1].push_back(i);
    }

    std::vector<std::vector<ConstraintBond>> fallback_bonds_by_color;
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
    std::vector<ConstraintBond> fallback_shake_bonds;
    for (int color = 0; color < fallback_bonds_by_color.size(); color++) {
        fallback_color_offsets[color] = fallback_shake_bonds.size();
        fallback_shake_bonds.insert(
            fallback_shake_bonds.end(),
            fallback_bonds_by_color[color].begin(),
            fallback_bonds_by_color[color].end());
    }
    fallback_color_offsets[fallback_bonds_by_color.size()] = fallback_shake_bonds.size();

    this->fallback_color_offsets = HostDeviceBuffer<int>::from_vector(fallback_color_offsets, true);
    fallback_n_colors = static_cast<int>(fallback_color_offsets.size()) - 1;
    this->fallback_shake_bonds = HostDeviceBuffer<ConstraintBond>::from_vector(fallback_shake_bonds, true);

    fallback_coop_blocks = 0;
    for (int c = 0; c < fallback_n_colors; c++) {
        const int n = fallback_color_offsets[c + 1] - fallback_color_offsets[c];
        fallback_coop_blocks = std::max(fallback_coop_blocks, (n + kShakeThreads - 1) / kShakeThreads);
    }
}

void CudaShake::apply(Context& ctx, HostDeviceBuffer<coord_t>& xcoords_buffer) {
    auto* coords = ctx.coords->gpu_data_p;
    auto* xcoords = xcoords_buffer.gpu_data_p;
    apply_to(ctx, coords, xcoords);
}

void CudaShake::initial_constraint(Context& ctx) {
    HostDeviceBuffer<coord_t> xcoords(ctx.n_atoms, true, true);
    auto* d_coords = ctx.coords->gpu_data_p;
    auto* d_xcoords = xcoords.gpu_data_p;

    check_cuda(cudaMemcpy(d_xcoords, d_coords, ctx.n_atoms * sizeof(coord_t), cudaMemcpyDeviceToDevice));
    apply_to(ctx, d_coords, d_xcoords);

    ctx.coords->download();

    auto* coords = ctx.coords->cpu_data_p;
    auto* velocities = ctx.velocities->cpu_data_p;
    auto* host_xcoords = xcoords.cpu_data_p;
    for (int i = 0; i < ctx.n_atoms; i++) {
        const double dt = ctx.dt;
        host_xcoords[i].x = coords[i].x - dt * velocities[i].x;
        host_xcoords[i].y = coords[i].y - dt * velocities[i].y;
        host_xcoords[i].z = coords[i].z - dt * velocities[i].z;
    }

    xcoords.upload();

    apply_to(ctx, d_xcoords, d_coords);

    xcoords.download();

    for (int i = 0; i < ctx.n_atoms; i++) {
        const double dt = ctx.dt;
        velocities[i].x = (coords[i].x - host_xcoords[i].x) / dt;
        velocities[i].y = (coords[i].y - host_xcoords[i].y) / dt;
        velocities[i].z = (coords[i].z - host_xcoords[i].z) / dt;
    }
    ctx.velocities->upload();
}

void CudaShake::init_backend(Context& ctx) {
    if (is_init_backend) {
        return;
    }

    std::vector<bool> optimized(data_.n_constraints, false);

    if (serial_q_molecules_) find_serial_q_molecule_bonds(ctx, optimized);
    find_shake_fast_water(ctx, optimized);
    find_shake_network(ctx, optimized);
    find_fallback_shake_bond(ctx, optimized);
    fallback_unconverged = std::make_unique<HostDeviceBuffer<int>>(1, true, true);
    shake_network_failed = std::make_unique<HostDeviceBuffer<int>>(1, true, true);
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

    const double xo0 = xcoords[o].x;
    const double yo0 = xcoords[o].y;
    const double zo0 = xcoords[o].z;
    const double xh10 = xcoords[h1].x;
    const double yh10 = xcoords[h1].y;
    const double zh10 = xcoords[h1].z;
    const double xh20 = xcoords[h2].x;
    const double yh20 = xcoords[h2].y;
    const double zh20 = xcoords[h2].z;

    const double xo1 = coords[o].x;
    const double yo1 = coords[o].y;
    const double zo1 = coords[o].z;
    const double xh11 = coords[h1].x;
    const double yh11 = coords[h1].y;
    const double zh11 = coords[h1].z;
    const double xh21 = coords[h2].x;
    const double yh21 = coords[h2].y;
    const double zh21 = coords[h2].z;

    const double xb0 = xh10 - xo0;
    const double yb0 = yh10 - yo0;
    const double zb0 = zh10 - zo0;
    const double xc0 = xh20 - xo0;
    const double yc0 = yh20 - yo0;
    const double zc0 = zh20 - zo0;

    const double xcom = xo1 * water.wo_div_wohh + (xh11 + xh21) * water.wh_div_wohh;
    const double ycom = yo1 * water.wo_div_wohh + (yh11 + yh21) * water.wh_div_wohh;
    const double zcom = zo1 * water.wo_div_wohh + (zh11 + zh21) * water.wh_div_wohh;

    const double xa1 = xo1 - xcom;
    const double ya1 = yo1 - ycom;
    const double za1 = zo1 - zcom;
    const double xb1 = xh11 - xcom;
    const double yb1 = yh11 - ycom;
    const double zb1 = zh11 - zcom;
    const double xc1 = xh21 - xcom;
    const double yc1 = yh21 - ycom;
    const double zc1 = zh21 - zcom;

    const double xakszd = yb0 * zc0 - zb0 * yc0;
    const double yakszd = zb0 * xc0 - xb0 * zc0;
    const double zakszd = xb0 * yc0 - yb0 * xc0;
    const double xaksxd = ya1 * zakszd - za1 * yakszd;
    const double yaksxd = za1 * xakszd - xa1 * zakszd;
    const double zaksxd = xa1 * yakszd - ya1 * xakszd;
    const double xaksyd = yakszd * zaksxd - zakszd * yaksxd;
    const double yaksyd = zakszd * xaksxd - xakszd * zaksxd;
    const double zaksyd = xakszd * yaksxd - yakszd * xaksxd;

    const double ax2 = xaksxd * xaksxd + yaksxd * yaksxd + zaksxd * zaksxd;
    const double ay2 = xaksyd * xaksyd + yaksyd * yaksyd + zaksyd * zaksyd;
    const double az2 = xakszd * xakszd + yakszd * yakszd + zakszd * zakszd;
    if (ax2 <= 0.0 || ay2 <= 0.0 || az2 <= 0.0) return;

    const double axlng_inv = rsqrt(ax2);
    const double aylng_inv = rsqrt(ay2);
    const double azlng_inv = rsqrt(az2);

    const double trns11 = xaksxd * axlng_inv;
    const double trns21 = yaksxd * axlng_inv;
    const double trns31 = zaksxd * axlng_inv;
    const double trns12 = xaksyd * aylng_inv;
    const double trns22 = yaksyd * aylng_inv;
    const double trns32 = zaksyd * aylng_inv;
    const double trns13 = xakszd * azlng_inv;
    const double trns23 = yakszd * azlng_inv;
    const double trns33 = zakszd * azlng_inv;

    const double xb0d = trns11 * xb0 + trns21 * yb0 + trns31 * zb0;
    const double yb0d = trns12 * xb0 + trns22 * yb0 + trns32 * zb0;
    const double xc0d = trns11 * xc0 + trns21 * yc0 + trns31 * zc0;
    const double yc0d = trns12 * xc0 + trns22 * yc0 + trns32 * zc0;
    const double za1d = trns13 * xa1 + trns23 * ya1 + trns33 * za1;
    const double xb1d = trns11 * xb1 + trns21 * yb1 + trns31 * zb1;
    const double yb1d = trns12 * xb1 + trns22 * yb1 + trns32 * zb1;
    const double zb1d = trns13 * xb1 + trns23 * yb1 + trns33 * zb1;
    const double xc1d = trns11 * xc1 + trns21 * yc1 + trns31 * zc1;
    const double yc1d = trns12 * xc1 + trns22 * yc1 + trns32 * zc1;
    const double zc1d = trns13 * xc1 + trns23 * yc1 + trns33 * zc1;

    const double sinphi = za1d * water.ra_inv;
    const double cosphi = sqrt(clamp_nonnegative(1.0 - sinphi * sinphi));
    if (cosphi <= 0.0) return;
    const double sinpsi = (zb1d - zc1d) / (water.rhh * cosphi);
    const double cospsi = sqrt(clamp_nonnegative(1.0 - sinpsi * sinpsi));

    const double ya2d = water.ra * cosphi;
    double xb2d = -water.rc * cospsi;
    const double yb2d = -water.rb * cosphi - water.rc * sinpsi * sinphi;
    const double yc2d = -water.rb * cosphi + water.rc * sinpsi * sinphi;
    xb2d = -0.5 *
           sqrt(clamp_nonnegative(water.rhh2 - (yb2d - yc2d) * (yb2d - yc2d) - (zb1d - zc1d) * (zb1d - zc1d)));

    const double alpa = xb2d * (xb0d - xc0d) + yb0d * yb2d + yc0d * yc2d;
    const double beta = xb2d * (yc0d - yb0d) + xb0d * yb2d + xc0d * yc2d;
    const double gama = xb0d * yb1d - xb1d * yb0d + xc0d * yc1d - xc1d * yc0d;
    const double al2be2 = alpa * alpa + beta * beta;
    if (al2be2 <= 0.0) return;

    const double sinthe = (alpa * gama - beta * sqrt(clamp_nonnegative(al2be2 - gama * gama))) / al2be2;
    const double costhe = sqrt(clamp_nonnegative(1.0 - sinthe * sinthe));

    const double xa3d = -ya2d * sinthe;
    const double ya3d = ya2d * costhe;
    const double za3d = za1d;
    const double xb3d = xb2d * costhe - yb2d * sinthe;
    const double yb3d = xb2d * sinthe + yb2d * costhe;
    const double zb3d = zb1d;
    const double xc3d = -xb2d * costhe - yc2d * sinthe;
    const double yc3d = -xb2d * sinthe + yc2d * costhe;
    const double zc3d = zc1d;

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
    coord_t* xcoords, int* failed) {
    int idx = blockIdx.x * blockDim.x + threadIdx.x;
    if (idx >= n_shake_networks) return;

    const auto network = networks[idx];
    double center_x = coords[network.center].x;
    double center_y = coords[network.center].y;
    double center_z = coords[network.center].z;
    double hydrogens_x[3];
    double hydrogens_y[3];
    double hydrogens_z[3];
    double old_vectors_x[3];
    double old_vectors_y[3];
    double old_vectors_z[3];

    for (int i = 0; i < network.n_hydrogens; i++) {
        const int h = network.hydrogens[i];
        hydrogens_x[i] = coords[h].x;
        hydrogens_y[i] = coords[h].y;
        hydrogens_z[i] = coords[h].z;
        old_vectors_x[i] = xcoords[network.center].x - xcoords[h].x;
        old_vectors_y[i] = xcoords[network.center].y - xcoords[h].y;
        old_vectors_z[i] = xcoords[network.center].z - xcoords[h].z;
    }

    bool converged = false;
    int n_iterations = 0;
    do {
        converged = true;
        for (int i = 0; i < network.n_hydrogens; i++) {
            const double current_vector_x = center_x - hydrogens_x[i];
            const double current_vector_y = center_y - hydrogens_y[i];
            const double current_vector_z = center_z - hydrogens_z[i];
            const double current_dist2 = current_vector_x * current_vector_x +
                                         current_vector_y * current_vector_y +
                                         current_vector_z * current_vector_z;
            const double target_dist2 = network.dist2[i];
            const double diff = target_dist2 - current_dist2;

            if (!isfinite(current_dist2) || !isfinite(diff)) {
                converged = false;
                continue;
            }

            if (fabs(diff) < shake_tol * network.dist2[i]) continue;

            converged = false;
            const double scp = current_vector_x * old_vectors_x[i] +
                               current_vector_y * old_vectors_y[i] +
                               current_vector_z * old_vectors_z[i];
            const double inv_mass_sum = network.center_winv + network.hydrogen_winv[i];
            if (!isfinite(scp) || !isfinite(inv_mass_sum) || scp <= kMinShakeProjectionFactor * target_dist2 || inv_mass_sum == 0.0) continue;

            const double corr = diff / (2.0 * scp * inv_mass_sum);
            const double center_scale = corr * network.center_winv;
            const double hydrogen_scale = corr * network.hydrogen_winv[i];

            if (!isfinite(corr) || !isfinite(center_scale) || !isfinite(hydrogen_scale) || fabs(center_scale) > kMaxShakeCorrectionScale || fabs(hydrogen_scale) > kMaxShakeCorrectionScale) {
                continue;
            }

            center_x += old_vectors_x[i] * center_scale;
            center_y += old_vectors_y[i] * center_scale;
            center_z += old_vectors_z[i] * center_scale;
            hydrogens_x[i] -= old_vectors_x[i] * hydrogen_scale;
            hydrogens_y[i] -= old_vectors_y[i] * hydrogen_scale;
            hydrogens_z[i] -= old_vectors_z[i] * hydrogen_scale;
        }
        n_iterations++;
    } while (n_iterations < shake_max_iter && !converged);

    if (!converged) {
        atomicExch(failed, 1);
        for (int i = 0; i < network.n_hydrogens; i++) {
            const double dx = center_x - hydrogens_x[i];
            const double dy = center_y - hydrogens_y[i];
            const double dz = center_z - hydrogens_z[i];
            const double dist2 = dx * dx + dy * dy + dz * dz;
            const double target_dist2 = network.dist2[i];

            if (!isfinite(dist2) || fabs(target_dist2 - dist2) >= shake_tol * target_dist2) {
                printf(">>> Shake failed, i = %d,j = %d, d = %f, d0 = %f\n",
                       network.center,
                       network.hydrogens[i],
                       sqrt(dist2),
                       sqrt(network.dist2[i]));
            }
        }
        return;
    }

    coords[network.center].x = center_x;
    coords[network.center].y = center_y;
    coords[network.center].z = center_z;
    for (int i = 0; i < network.n_hydrogens; i++) {
        coords[network.hydrogens[i]].x = hydrogens_x[i];
        coords[network.hydrogens[i]].y = hydrogens_y[i];
        coords[network.hydrogens[i]].z = hydrogens_z[i];
    }
}

__global__ void print_fallback_shake_failures_kernel(
    int n_shakes,
    ConstraintBond* shake_bonds,
    coord_t* coords) {
    const int idx = blockIdx.x * blockDim.x + threadIdx.x;
    if (idx >= n_shakes) return;

    const int ai = shake_bonds[idx].ai - 1;
    const int aj = shake_bonds[idx].aj - 1;
    const double dx = coords[ai].x - coords[aj].x;
    const double dy = coords[ai].y - coords[aj].y;
    const double dz = coords[ai].z - coords[aj].z;
    const double dist2 = dx * dx + dy * dy + dz * dz;
    if (!isfinite(dist2) || fabs(shake_bonds[idx].dist2 - dist2) >= shake_tol * shake_bonds[idx].dist2) {
        printf(">>> Shake failed, i = %d,j = %d, d = %f, d0 = %f\n",
               ai,
               aj,
               sqrt(dist2),
               sqrt(shake_bonds[idx].dist2));
    }
}

__global__ void fallback_shake_fused_kernel(
    int n_colors,
    const int* color_offsets,
    ConstraintBond* shake_bonds,
    coord_t* coords,
    coord_t* xcoords,
    double* winv,
    int* unconverged) {
    cg::grid_group grid = cg::this_grid();
    const int tid = blockIdx.x * blockDim.x + threadIdx.x;

    for (int iter = 0; iter < shake_max_iter; iter++) {
        if (grid.thread_rank() == 0) *unconverged = 0;
        grid.sync();

        for (int c = 0; c < n_colors; c++) {
            const int off = color_offsets[c];
            const int n = color_offsets[c + 1] - off;

            if (tid < n) {
                ConstraintBond& shake_bond = shake_bonds[off + tid];
                const int ai = shake_bond.ai - 1;
                const int aj = shake_bond.aj - 1;
                const double xij_x = coords[ai].x - coords[aj].x;
                const double xij_y = coords[ai].y - coords[aj].y;
                const double xij_z = coords[ai].z - coords[aj].z;
                const double xij2 = xij_x * xij_x + xij_y * xij_y + xij_z * xij_z;
                const double diff = shake_bond.dist2 - xij2;
                const double target_dist2 = shake_bond.dist2;
                const double tolerance = shake_tol * target_dist2;
                if (!isfinite(xij2) || !isfinite(diff)) {
                    *unconverged = 1;
                }

                if (fabs(diff) >= tolerance) {
                    *unconverged = 1;  // plain store is fine (any thread -> 1)
                    const double xxij_x = xcoords[ai].x - xcoords[aj].x;
                    const double xxij_y = xcoords[ai].y - xcoords[aj].y;
                    const double xxij_z = xcoords[ai].z - xcoords[aj].z;
                    const double scp = xij_x * xxij_x + xij_y * xxij_y + xij_z * xxij_z;
                    const double inv_mass_sum = winv[ai] + winv[aj];
                    const double min_scp = kMinShakeProjectionFactor * target_dist2;
                    const bool valid_projection =
                        isfinite(scp) &&
                        isfinite(inv_mass_sum) &&
                        scp > min_scp &&
                        inv_mass_sum > 0.0;
                    if (valid_projection) {
                        const double corr = diff / (2.0 * scp * inv_mass_sum);
                        const double ai_scale = corr * winv[ai];
                        const double aj_scale = corr * winv[aj];

                        const bool valid_correction = isfinite(corr) &&
                                                      isfinite(ai_scale) &&
                                                      isfinite(aj_scale) &&
                                                      fabs(ai_scale) <= kMaxShakeCorrectionScale &&
                                                      fabs(aj_scale) <= kMaxShakeCorrectionScale;

                        if (valid_correction) {
                            coords[ai].x += xxij_x * ai_scale;
                            coords[ai].y += xxij_y * ai_scale;
                            coords[ai].z += xxij_z * ai_scale;
                            coords[aj].x -= xxij_x * aj_scale;
                            coords[aj].y -= xxij_y * aj_scale;
                            coords[aj].z -= xxij_z * aj_scale;
                        }
                    }
                }
            }
            grid.sync();
        }
        if (*unconverged == 0) break;
        grid.sync();
    }
}

void CudaShake::apply_to(Context& ctx, coord_t* d_coords, coord_t* d_xcoords) {
    serial_q_solver_.apply(d_coords, d_xcoords, ctx.winv->gpu_data_p);
    if (shake_fast_waters->length > 0) {
        const int grid_blocks = (shake_fast_waters->length + kShakeThreads - 1) / kShakeThreads;
        calc_fast_water_shake_kernel<<<grid_blocks, kShakeThreads>>>(shake_fast_waters->length, shake_fast_waters->gpu_data_p, d_coords, d_xcoords);
    }
    if (shake_networks->length > 0) {
        const int grid_blocks = (shake_networks->length + kShakeThreads - 1) / kShakeThreads;
        shake_network_failed->zero();
        calc_h_star_shake_kernel<<<grid_blocks, kShakeThreads>>>(shake_networks->length, shake_networks->gpu_data_p, d_coords, d_xcoords, shake_network_failed->gpu_data_p);
        check_cuda(cudaGetLastError());
        shake_network_failed->download();

        if (shake_network_failed->cpu_data_p[0]) {
            std::fflush(stdout);
            std::exit(EXIT_FAILURE);
        }
    }
    if (fallback_shake_bonds->length > 0) {
        fallback_unconverged->zero();
        void* args[] = {
            &fallback_n_colors,
            &fallback_color_offsets->gpu_data_p,
            &fallback_shake_bonds->gpu_data_p,
            &d_coords,
            &d_xcoords,
            &ctx.winv->gpu_data_p,
            &fallback_unconverged->gpu_data_p,
        };
        dim3 grid(fallback_coop_blocks), block(kShakeThreads);
        check_cuda(cudaLaunchCooperativeKernel((void*)fallback_shake_fused_kernel, grid, block, args));
        fallback_unconverged->download();
        if (fallback_unconverged->cpu_data_p[0]) {
            int n_fallback_constraints = fallback_shake_bonds->length;
            const int grid_blocks = (n_fallback_constraints + kShakeThreads - 1) / kShakeThreads;
            print_fallback_shake_failures_kernel<<<grid_blocks, kShakeThreads>>>(
                n_fallback_constraints,
                fallback_shake_bonds->gpu_data_p,
                d_coords);
            check_cuda(cudaGetLastError());
            check_cuda(cudaDeviceSynchronize());
            std::fflush(stdout);
            std::exit(EXIT_FAILURE);
        }
    }
}

void CudaShake::cleanup() {
}
