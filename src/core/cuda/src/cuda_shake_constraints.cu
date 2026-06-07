#include <cooperative_groups.h>

#include <algorithm>
#include <cmath>
#include <cstdlib>
#include <iostream>
#include <unordered_map>
#include <vector>

#include "common/include/constants.h"
#include "common/include/context.h"
#include "cuda/include/cuda_shake_constraints.cuh"

namespace CudaShakeConstraints {

bool is_initialized = false;
constexpr int kShakeThreads = 64;
constexpr int kFallbackThreads = 256;
constexpr int kFallbackCooperativeThreads = 64;

struct shake_fast_water_t {
    int o;
    int h1;
    int h2;
    real_t ra;
    real_t ra_inv;
    real_t rb;
    real_t rc;
    real_t rc2;
    real_t hhhh;
    real_t wo_div_wohh;
    real_t wh_div_wohh;
};

struct shake_network_t {
    int center;
    int n_hydrogens;
    int hydrogens[3];
    real_t dist2[3];
    real_t center_winv;
    real_t hydrogen_winv[3];
};

int n_fast_waters = 0;
int n_shake_networks = 0;
int n_fallback_constraints = 0;
int n_fallback_components = 0;
int n_fallback_colors = 0;
shake_fast_water_t* d_fast_waters = nullptr;
shake_network_t* d_shake_networks = nullptr;
ShakeBond* d_fallback_shake_bonds = nullptr;
int* d_fallback_unconverged = nullptr;
int* d_fallback_color_offsets = nullptr;
cudaGraphExec_t fallback_sweep_graph_exec = nullptr;
coord_t* fallback_graph_coords = nullptr;
coord_t* fallback_graph_xcoords = nullptr;
real_t* fallback_graph_winv = nullptr;
bool fallback_launch_config_initialized = false;
bool fallback_use_cooperative = false;
int fallback_cooperative_grid_blocks = 0;
std::vector<int> h_fallback_color_offsets;

long long pair_key(int a, int b) {
    if (a > b) std::swap(a, b);
    return (static_cast<long long>(a) << 32) | static_cast<unsigned int>(b);
}

bool nearly_equal(real_t a, real_t b) {
    const real_t scale = std::max(static_cast<real_t>(1), std::max(std::abs(a), std::abs(b)));
    return std::abs(a - b) <= static_cast<real_t>(1.0e-5) * scale;
}

int lookup_pair(const std::unordered_map<long long, int>& pair_to_bond, int a, int b) {
    auto it = pair_to_bond.find(pair_key(a, b));
    if (it == pair_to_bond.end()) return -1;
    return it->second;
}

shake_fast_water_t make_fast_water(int o, int h1, int h2, real_t roh2, real_t rhh2, real_t wo, real_t wh) {
    shake_fast_water_t water = {};
    water.o = o;
    water.h1 = h1;
    water.h2 = h2;

    const real_t roh = std::sqrt(roh2);
    const real_t rhh = std::sqrt(rhh2);
    const real_t height = std::sqrt(std::max(static_cast<real_t>(0), roh * roh - static_cast<real_t>(0.25) * rhh * rhh));
    const real_t total_mass = wo + static_cast<real_t>(2) * wh;
    const real_t com_y = -static_cast<real_t>(2) * height * wh / total_mass;
    water.ra = std::abs(com_y);
    water.ra_inv = static_cast<real_t>(1) / water.ra;
    water.rb = height - water.ra;
    water.rc = static_cast<real_t>(0.5) * rhh;
    water.rc2 = rhh;
    water.hhhh = rhh * rhh;
    water.wo_div_wohh = wo / total_mass;
    water.wh_div_wohh = wh / total_mass;
    return water;
}
}  // namespace CudaShakeConstraints

__device__ real_t clamp_nonnegative(real_t value) {
    return value > static_cast<real_t>(0) ? value : static_cast<real_t>(0);
}

__global__ void calc_fast_water_shake_kernel(
    int n_fast_waters,
    CudaShakeConstraints::shake_fast_water_t* fast_waters,
    coord_t* coords,
    coord_t* xcoords) {
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
    const real_t sinpsi = (zb1d - zc1d) / (water.rc2 * cosphi);
    const real_t cospsi = sqrt(clamp_nonnegative(static_cast<real_t>(1) - sinpsi * sinpsi));

    const real_t ya2d = water.ra * cosphi;
    real_t xb2d = -water.rc * cospsi;
    const real_t yb2d = -water.rb * cosphi - water.rc * sinpsi * sinphi;
    const real_t yc2d = -water.rb * cosphi + water.rc * sinpsi * sinphi;
    xb2d = -static_cast<real_t>(0.5) *
           sqrt(clamp_nonnegative(water.hhhh - (yb2d - yc2d) * (yb2d - yc2d) - (zb1d - zc1d) * (zb1d - zc1d)));

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
    CudaShakeConstraints::shake_network_t* networks,
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

__device__ void apply_fallback_shake_bond(
    const ShakeBond& shake_bond,
    coord_t* coords,
    coord_t* xcoords,
    real_t* winv,
    int* unconverged) {
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

__global__ void calc_fallback_shake_color_kernel(
    int n_shakes,
    ShakeBond* shake_bonds,
    coord_t* coords,
    coord_t* xcoords,
    real_t* winv,
    int* unconverged) {
    const int idx = blockIdx.x * blockDim.x + threadIdx.x;
    if (idx >= n_shakes) return;

    apply_fallback_shake_bond(shake_bonds[idx], coords, xcoords, winv, unconverged);
}

__global__ void calc_fallback_shake_cooperative_kernel(
    int n_colors,
    int* color_offsets,
    ShakeBond* shake_bonds,
    coord_t* coords,
    coord_t* xcoords,
    real_t* winv,
    int* unconverged) {
    auto grid = cooperative_groups::this_grid();
    const int thread_rank = static_cast<int>(grid.thread_rank());
    const int grid_threads = static_cast<int>(grid.size());
    bool converged = false;

    for (int n_iterations = 0; n_iterations < shake_max_iter; n_iterations++) {
        if (thread_rank == 0) *unconverged = 0;
        grid.sync();

        for (int color = 0; color < n_colors; color++) {
            const int offset = color_offsets[color];
            const int n_color_shakes = color_offsets[color + 1] - offset;
            for (int idx = thread_rank; idx < n_color_shakes; idx += grid_threads) {
                apply_fallback_shake_bond(shake_bonds[offset + idx], coords, xcoords, winv, unconverged);
            }
            grid.sync();
        }

        converged = (*unconverged == 0);
        grid.sync();
        if (converged) break;
    }

    if (!converged) {
        const int n_shakes = color_offsets[n_colors];
        for (int idx = thread_rank; idx < n_shakes; idx += grid_threads) {
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

void destroy_fallback_sweep_graph() {
    using namespace CudaShakeConstraints;
    if (fallback_sweep_graph_exec) {
        check_cuda(cudaGraphExecDestroy(fallback_sweep_graph_exec));
        fallback_sweep_graph_exec = nullptr;
    }
    fallback_graph_coords = nullptr;
    fallback_graph_xcoords = nullptr;
    fallback_graph_winv = nullptr;
}

void ensure_fallback_sweep_graph(coord_t* d_coords, coord_t* d_xcoords, real_t* d_winv) {
    using namespace CudaShakeConstraints;
    if (n_fallback_constraints == 0) return;
    if (fallback_sweep_graph_exec &&
        fallback_graph_coords == d_coords &&
        fallback_graph_xcoords == d_xcoords &&
        fallback_graph_winv == d_winv) {
        return;
    }

    destroy_fallback_sweep_graph();

    cudaStream_t capture_stream = nullptr;
    cudaGraph_t graph = nullptr;
    check_cuda(cudaStreamCreate(&capture_stream));
    check_cuda(cudaStreamBeginCapture(capture_stream, cudaStreamCaptureModeGlobal));

    check_cuda(cudaMemsetAsync(d_fallback_unconverged, 0, sizeof(int), capture_stream));
    for (int color = 0; color < n_fallback_colors; color++) {
        const int offset = h_fallback_color_offsets[color];
        const int n_color_shakes = h_fallback_color_offsets[color + 1] - offset;
        if (n_color_shakes == 0) continue;
        const int grid_blocks = (n_color_shakes + kFallbackThreads - 1) / kFallbackThreads;
        calc_fallback_shake_color_kernel<<<grid_blocks, kFallbackThreads, 0, capture_stream>>>(
            n_color_shakes,
            d_fallback_shake_bonds + offset,
            d_coords,
            d_xcoords,
            d_winv,
            d_fallback_unconverged);
    }

    check_cuda(cudaStreamEndCapture(capture_stream, &graph));
    check_cuda(cudaGraphInstantiate(&fallback_sweep_graph_exec, graph, nullptr, nullptr, 0));
    check_cuda(cudaGraphDestroy(graph));
    check_cuda(cudaStreamDestroy(capture_stream));

    fallback_graph_coords = d_coords;
    fallback_graph_xcoords = d_xcoords;
    fallback_graph_winv = d_winv;
}

void ensure_fallback_launch_config() {
    using namespace CudaShakeConstraints;
    if (fallback_launch_config_initialized) return;
    fallback_launch_config_initialized = true;
    fallback_use_cooperative = false;
    fallback_cooperative_grid_blocks = 0;

    if (n_fallback_constraints == 0 || n_fallback_colors == 0 || d_fallback_color_offsets == nullptr) return;

    int device = 0;
    int cooperative_launch = 0;
    int sm_count = 0;
    int max_active_blocks_per_sm = 0;
    check_cuda(cudaGetDevice(&device));
    check_cuda(cudaDeviceGetAttribute(&cooperative_launch, cudaDevAttrCooperativeLaunch, device));
    check_cuda(cudaDeviceGetAttribute(&sm_count, cudaDevAttrMultiProcessorCount, device));
    if (!cooperative_launch || sm_count <= 0) return;

    check_cuda(cudaOccupancyMaxActiveBlocksPerMultiprocessor(
        &max_active_blocks_per_sm,
        calc_fallback_shake_cooperative_kernel,
        kFallbackCooperativeThreads,
        0));
    if (max_active_blocks_per_sm <= 0) return;

    int max_color_shakes = 0;
    for (int color = 0; color < n_fallback_colors; color++) {
        const int n_color_shakes = h_fallback_color_offsets[color + 1] - h_fallback_color_offsets[color];
        max_color_shakes = std::max(max_color_shakes, n_color_shakes);
    }
    if (max_color_shakes <= 0) return;

    const int requested_blocks = (max_color_shakes + kFallbackCooperativeThreads - 1) / kFallbackCooperativeThreads;
    const int max_cooperative_blocks = sm_count * max_active_blocks_per_sm;
    fallback_cooperative_grid_blocks = std::max(1, std::min(requested_blocks, max_cooperative_blocks));
    fallback_use_cooperative = fallback_cooperative_grid_blocks > 0;

    if (std::getenv("QDYN_SHAKE_DEBUG") != nullptr) {
        std::cerr << "CUDA SHAKE fallback launch: mode="
                  << (fallback_use_cooperative ? "cooperative" : "graph")
                  << ", cooperative_grid_blocks=" << fallback_cooperative_grid_blocks
                  << ", cooperative_threads=" << kFallbackCooperativeThreads
                  << ", max_color_shakes=" << max_color_shakes
                  << ", max_cooperative_blocks=" << max_cooperative_blocks << std::endl;
    }
}

bool launch_fallback_shake_cooperative(coord_t* d_coords, coord_t* d_xcoords, real_t* d_winv) {
    using namespace CudaShakeConstraints;
    ensure_fallback_launch_config();
    if (!fallback_use_cooperative) return false;

    void* args[] = {
        &n_fallback_colors,
        &d_fallback_color_offsets,
        &d_fallback_shake_bonds,
        &d_coords,
        &d_xcoords,
        &d_winv,
        &d_fallback_unconverged,
    };
    check_cuda(cudaLaunchCooperativeKernel(
        reinterpret_cast<void*>(calc_fallback_shake_cooperative_kernel),
        fallback_cooperative_grid_blocks,
        kFallbackCooperativeThreads,
        args,
        0,
        0));
    return true;
}

void init_shake_constraints_kernel_data(Context& ctx) {
    auto& host = ctx;
    using namespace CudaShakeConstraints;
    if (!is_initialized) {
        if (host.shake_data.n_constraints > 0) {
            auto* shake_bonds = host.shake_data.shake_bonds->cpu_data_p;
            auto* heavy = host.heavy->cpu_data_p;
            auto* winv = host.winv->cpu_data_p;

            std::vector<int> bond_ai(host.shake_data.n_constraints);
            std::vector<int> bond_aj(host.shake_data.n_constraints);
            std::vector<int> atom_degree(host.n_atoms, 0);
            std::unordered_map<long long, int> pair_to_bond;
            pair_to_bond.reserve(static_cast<size_t>(host.shake_data.n_constraints) * 2);

            for (int i = 0; i < host.shake_data.n_constraints; i++) {
                const int ai = shake_bonds[i].ai - 1;
                const int aj = shake_bonds[i].aj - 1;
                bond_ai[i] = ai;
                bond_aj[i] = aj;
                if (ai >= 0 && ai < host.n_atoms) atom_degree[ai]++;
                if (aj >= 0 && aj < host.n_atoms) atom_degree[aj]++;
                if (ai >= 0 && ai < host.n_atoms && aj >= 0 && aj < host.n_atoms) {
                    const long long key = pair_key(ai, aj);
                    auto [it, inserted] = pair_to_bond.emplace(key, i);
                    if (!inserted) it->second = -1;
                }
            }

            std::vector<bool> optimized(host.shake_data.n_constraints, false);
            std::vector<shake_fast_water_t> fast_waters;
            fast_waters.reserve(host.n_waters);

            for (int w = 0; w < host.n_waters; w++) {
                const int o = host.n_atoms_solute + 3 * w;
                const int h1 = o + 1;
                const int h2 = o + 2;
                if (h2 >= host.n_atoms) continue;
                if (!heavy[o] || heavy[h1] || heavy[h2]) continue;
                if (atom_degree[o] != 2 || atom_degree[h1] != 2 || atom_degree[h2] != 2) continue;

                const int oh1 = lookup_pair(pair_to_bond, o, h1);
                const int oh2 = lookup_pair(pair_to_bond, o, h2);
                const int hh = lookup_pair(pair_to_bond, h1, h2);
                if (oh1 < 0 || oh2 < 0 || hh < 0) continue;
                if (!nearly_equal(shake_bonds[oh1].dist2, shake_bonds[oh2].dist2)) continue;
                if (shake_bonds[oh1].dist2 <= static_cast<real_t>(0) || shake_bonds[hh].dist2 <= static_cast<real_t>(0)) continue;
                if (shake_bonds[oh1].dist2 <= static_cast<real_t>(0.25) * shake_bonds[hh].dist2) continue;
                if (winv[o] <= static_cast<real_t>(0) || winv[h1] <= static_cast<real_t>(0) || winv[h2] <= static_cast<real_t>(0)) continue;

                const real_t wo = static_cast<real_t>(1) / winv[o];
                const real_t wh1 = static_cast<real_t>(1) / winv[h1];
                const real_t wh2 = static_cast<real_t>(1) / winv[h2];
                if (!nearly_equal(wh1, wh2)) continue;

                fast_waters.push_back(make_fast_water(o, h1, h2, shake_bonds[oh1].dist2, shake_bonds[hh].dist2, wo, wh1));
                optimized[oh1] = true;
                optimized[oh2] = true;
                optimized[hh] = true;
            }

            std::vector<std::vector<int>> atom_to_bonds(host.n_atoms);
            for (int i = 0; i < host.shake_data.n_constraints; i++) {
                if (optimized[i]) continue;
                if (bond_ai[i] >= 0 && bond_ai[i] < host.n_atoms) atom_to_bonds[bond_ai[i]].push_back(i);
                if (bond_aj[i] >= 0 && bond_aj[i] < host.n_atoms) atom_to_bonds[bond_aj[i]].push_back(i);
            }

            std::vector<shake_network_t> shake_networks;
            std::vector<bool> visited(host.shake_data.n_constraints, false);
            std::vector<int> atom_seen(host.n_atoms, -1);
            int component_id = 0;

            for (int start = 0; start < host.shake_data.n_constraints; start++) {
                if (optimized[start] || visited[start]) continue;

                std::vector<int> component_bonds;
                std::vector<int> component_atoms;
                std::vector<int> stack = {start};
                visited[start] = true;

                while (!stack.empty()) {
                    const int bi = stack.back();
                    stack.pop_back();
                    component_bonds.push_back(bi);

                    const int atoms[2] = {bond_ai[bi], bond_aj[bi]};
                    for (int atom : atoms) {
                        if (atom < 0 || atom >= host.n_atoms) continue;
                        if (atom_seen[atom] != component_id) {
                            atom_seen[atom] = component_id;
                            component_atoms.push_back(atom);
                        }
                        for (int next_bond : atom_to_bonds[atom]) {
                            if (!visited[next_bond]) {
                                visited[next_bond] = true;
                                stack.push_back(next_bond);
                            }
                        }
                    }
                }
                component_id++;

                if (component_bonds.empty() || component_bonds.size() > 3 || component_atoms.size() != component_bonds.size() + 1) continue;

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

                shake_network_t network = {};
                network.center = center;
                network.center_winv = winv[center];
                bool valid = true;
                for (int bi : component_bonds) {
                    const int ai = bond_ai[bi];
                    const int aj = bond_aj[bi];
                    int hydrogen = -1;
                    if (ai == center && aj >= 0 && aj < host.n_atoms && !heavy[aj]) {
                        hydrogen = aj;
                    } else if (aj == center && ai >= 0 && ai < host.n_atoms && !heavy[ai]) {
                        hydrogen = ai;
                    } else {
                        valid = false;
                        break;
                    }

                    const int hidx = network.n_hydrogens;
                    if (hidx >= 3 || network.center_winv + winv[hydrogen] == static_cast<real_t>(0)) {
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

            std::vector<std::vector<int>> fallback_atom_to_bonds(host.n_atoms);
            for (int i = 0; i < host.shake_data.n_constraints; i++) {
                if (optimized[i]) continue;
                if (bond_ai[i] >= 0 && bond_ai[i] < host.n_atoms) fallback_atom_to_bonds[bond_ai[i]].push_back(i);
                if (bond_aj[i] >= 0 && bond_aj[i] < host.n_atoms) fallback_atom_to_bonds[bond_aj[i]].push_back(i);
            }

            std::vector<std::vector<ShakeBond>> fallback_bonds_by_color;
            std::vector<ShakeBond> fallback_shake_bonds;
            fallback_shake_bonds.reserve(host.shake_data.n_constraints);
            int max_fallback_colors = 0;
            int fallback_component_count = 0;
            std::vector<bool> fallback_visited(host.shake_data.n_constraints, false);
            for (int start = 0; start < host.shake_data.n_constraints; start++) {
                if (optimized[start] || fallback_visited[start]) continue;

                std::vector<int> component_bonds;
                std::vector<int> stack = {start};
                fallback_visited[start] = true;
                while (!stack.empty()) {
                    const int bi = stack.back();
                    stack.pop_back();
                    component_bonds.push_back(bi);

                    const int atoms[2] = {bond_ai[bi], bond_aj[bi]};
                    for (int atom : atoms) {
                        if (atom < 0 || atom >= host.n_atoms) continue;
                        for (int next_bond : fallback_atom_to_bonds[atom]) {
                            if (!fallback_visited[next_bond]) {
                                fallback_visited[next_bond] = true;
                                stack.push_back(next_bond);
                            }
                        }
                    }
                }

                std::unordered_map<int, std::vector<int>> atom_to_colors;
                atom_to_colors.reserve(component_bonds.size() * 2);
                std::vector<int> component_colors;
                component_colors.reserve(component_bonds.size());
                int component_n_colors = 0;
                for (int bi : component_bonds) {
                    std::vector<char> used(component_n_colors + 1, 0);
                    const int ai = bond_ai[bi];
                    const int aj = bond_aj[bi];
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

                fallback_component_count++;
                max_fallback_colors = std::max(max_fallback_colors, component_n_colors);
                if (component_n_colors > static_cast<int>(fallback_bonds_by_color.size())) {
                    fallback_bonds_by_color.resize(component_n_colors);
                }
                for (size_t i = 0; i < component_bonds.size(); i++) {
                    const int bi = component_bonds[i];
                    fallback_bonds_by_color[component_colors[i]].push_back(shake_bonds[bi]);
                }
            }

            std::vector<int> fallback_color_offsets(fallback_bonds_by_color.size() + 1, 0);
            for (size_t color = 0; color < fallback_bonds_by_color.size(); color++) {
                fallback_color_offsets[color] = static_cast<int>(fallback_shake_bonds.size());
                fallback_shake_bonds.insert(
                    fallback_shake_bonds.end(),
                    fallback_bonds_by_color[color].begin(),
                    fallback_bonds_by_color[color].end());
            }
            fallback_color_offsets[fallback_bonds_by_color.size()] = static_cast<int>(fallback_shake_bonds.size());
            h_fallback_color_offsets = fallback_color_offsets;

            n_fast_waters = static_cast<int>(fast_waters.size());
            n_shake_networks = static_cast<int>(shake_networks.size());
            n_fallback_constraints = static_cast<int>(fallback_shake_bonds.size());
            n_fallback_components = fallback_component_count;
            n_fallback_colors = static_cast<int>(fallback_bonds_by_color.size());

            if (std::getenv("QDYN_SHAKE_DEBUG") != nullptr) {
                std::cerr << "CUDA SHAKE: fast_waters=" << n_fast_waters
                          << ", h_star_networks=" << n_shake_networks
                          << ", fallback_components=" << n_fallback_components
                          << ", fallback_constraints=" << n_fallback_constraints
                          << ", fallback_colors=" << n_fallback_colors
                          << ", max_component_colors=" << max_fallback_colors << std::endl;
            }

            if (n_fast_waters > 0) {
                check_cudaMalloc((void**)&d_fast_waters, sizeof(shake_fast_water_t) * n_fast_waters);
                check_cuda(cudaMemcpy(d_fast_waters, fast_waters.data(), sizeof(shake_fast_water_t) * n_fast_waters, cudaMemcpyHostToDevice));
            }
            if (n_shake_networks > 0) {
                check_cudaMalloc((void**)&d_shake_networks, sizeof(shake_network_t) * n_shake_networks);
                check_cuda(cudaMemcpy(d_shake_networks, shake_networks.data(), sizeof(shake_network_t) * n_shake_networks, cudaMemcpyHostToDevice));
            }
            if (n_fallback_constraints > 0) {
                check_cudaMalloc((void**)&d_fallback_shake_bonds, sizeof(ShakeBond) * n_fallback_constraints);
                check_cudaMalloc((void**)&d_fallback_unconverged, sizeof(int));
                check_cudaMalloc((void**)&d_fallback_color_offsets, sizeof(int) * h_fallback_color_offsets.size());
                check_cuda(cudaMemcpy(d_fallback_shake_bonds, fallback_shake_bonds.data(), sizeof(ShakeBond) * n_fallback_constraints, cudaMemcpyHostToDevice));
                check_cuda(cudaMemcpy(d_fallback_color_offsets, h_fallback_color_offsets.data(), sizeof(int) * h_fallback_color_offsets.size(), cudaMemcpyHostToDevice));
            }
        }

        is_initialized = true;
    }
}

void cleanup_shake_constraints() {
    using namespace CudaShakeConstraints;
    if (is_initialized) {
        destroy_fallback_sweep_graph();
        if (d_fast_waters) cudaFree(d_fast_waters);
        if (d_shake_networks) cudaFree(d_shake_networks);
        if (d_fallback_shake_bonds) cudaFree(d_fallback_shake_bonds);
        if (d_fallback_unconverged) cudaFree(d_fallback_unconverged);
        if (d_fallback_color_offsets) cudaFree(d_fallback_color_offsets);
        d_fast_waters = nullptr;
        d_shake_networks = nullptr;
        d_fallback_shake_bonds = nullptr;
        d_fallback_unconverged = nullptr;
        d_fallback_color_offsets = nullptr;
        n_fast_waters = 0;
        n_shake_networks = 0;
        n_fallback_constraints = 0;
        n_fallback_components = 0;
        n_fallback_colors = 0;
        fallback_launch_config_initialized = false;
        fallback_use_cooperative = false;
        fallback_cooperative_grid_blocks = 0;
        h_fallback_color_offsets.clear();
        is_initialized = false;
    }
}

void calc_shake_constraints_host(Context& ctx) {
    auto& host = ctx;
    if (host.n_molecules == 0 || host.shake_data.n_constraints == 0) return;
    using namespace CudaShakeConstraints;

    auto d_coords = host.coords->gpu_data_p;
    auto d_xcoords = host.xcoords->gpu_data_p;
    auto d_winv = host.winv->gpu_data_p;
    printf("xxxxxxxxxxxxxxxxxxxxxxxx %d %d %d\n", n_fast_waters, n_shake_networks, n_fallback_constraints);

    if (n_fast_waters > 0) {
        const int grid_blocks = (n_fast_waters + kShakeThreads - 1) / kShakeThreads;
        calc_fast_water_shake_kernel<<<grid_blocks, kShakeThreads>>>(n_fast_waters, d_fast_waters, d_coords, d_xcoords);
    }
    if (n_shake_networks > 0) {
        const int grid_blocks = (n_shake_networks + kShakeThreads - 1) / kShakeThreads;
        calc_h_star_shake_kernel<<<grid_blocks, kShakeThreads>>>(n_shake_networks, d_shake_networks, d_coords, d_xcoords);
    }
    if (n_fallback_constraints > 0) {
        if (launch_fallback_shake_cooperative(d_coords, d_xcoords, d_winv)) return;

        ensure_fallback_sweep_graph(d_coords, d_xcoords, d_winv);

        int unconverged_host = 0;
        bool converged = false;
        for (int n_iterations = 0; n_iterations < shake_max_iter; n_iterations++) {
            check_cuda(cudaGraphLaunch(fallback_sweep_graph_exec, 0));
            check_cuda(cudaMemcpy(&unconverged_host, d_fallback_unconverged, sizeof(int), cudaMemcpyDeviceToHost));
            if (unconverged_host == 0) {
                converged = true;
                break;
            }
        }

        if (!converged) {
            const int grid_blocks = (n_fallback_constraints + kFallbackThreads - 1) / kFallbackThreads;
            print_fallback_shake_failures_kernel<<<grid_blocks, kFallbackThreads>>>(
                n_fallback_constraints,
                d_fallback_shake_bonds,
                d_coords);
        }
    }
}
