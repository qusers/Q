#include <iostream>

#include "context.h"
#include "cuda_context.cuh"
#include "cuda_nonbonded_force.cuh"
#include "constants.h"
#include "vdw_rules.h"
#include "cuda_utility.cuh"

namespace CudaNonbondedForce {
bool is_initialized = false;
double *d_evdw_total, *d_ecoul_total;

__device__ __forceinline__ void idx2xy(int n, int t, int& x, int& y) {
    x = (int)floorf((2 * n + 1 - sqrtf((2 * n + 1) * (2 * n + 1) - 8 * t)) * 0.5f);
    y = t - (x * n - (x * (x - 1) >> 1));
    if (y < 0) {
        x--;
        y += (n - x);
    }
    y += x;
}

__device__ __forceinline__ double shfl(double v, int srcLane, unsigned mask = 0xffffffffu) {
    int2 a = *reinterpret_cast<int2*>(&v);
    a.x = __shfl_sync(mask, a.x, srcLane);
    a.y = __shfl_sync(mask, a.y, srcLane);
    return *reinterpret_cast<double*>(&a);
}

__device__ __forceinline__ coord_t shfl_coord(coord_t v, int srcLane, unsigned mask = 0xffffffffu) {
    v.x = shfl(v.x, srcLane, mask);
    v.y = shfl(v.y, srcLane, mask);
    v.z = shfl(v.z, srcLane, mask);
    return v;
}

__device__ __forceinline__ catype_t shfl_catype(catype_t v, int srcLane, unsigned mask = 0xffffffffu) {
    v.code = __shfl_sync(mask, v.code, srcLane);
    v.m = shfl(v.m, srcLane, mask);
    v.aii_normal = shfl(v.aii_normal, srcLane, mask);
    v.bii_normal = shfl(v.bii_normal, srcLane, mask);
    v.aii_polar = shfl(v.aii_polar, srcLane, mask);
    v.bii_polar = shfl(v.bii_polar, srcLane, mask);
    v.aii_1_4 = shfl(v.aii_1_4, srcLane, mask);
    v.bii_1_4 = shfl(v.bii_1_4, srcLane, mask);
    return v;
}

__device__ void calculate_unforce_bound(
    const coord_t& x,
    const coord_t& y,

    const double x_charge,
    const double y_charge,

    const double x_aii,
    const double y_aii,

    const double x_bii,
    const double y_bii,

    const double coulomb_constant,

    const double scaling,

    const int vdw_rule,
    const double lambda,

    double& evdw,
    double& ecoul,
    double& dv) {
    double3 d = {x.x - y.x, x.y - y.y, x.z - y.z};
    double r = rsqrt(d.x * d.x + d.y * d.y + d.z * d.z);
    double r2 = r * r;
    double r6 = r2 * r2 * r2;
    // double v_a = r6 * r6;
    // double v_b = r6;
    // ecoul = r;
    // evdw = v_a - v_b;
    // dv = r2 * (-ecoul - v_a + v_b);

    ecoul = scaling * coulomb_constant * x_charge * y_charge * r * lambda;

    double v_a, v_b;
    if (vdw_rule == VDW_GEOMETRIC) {
        calc_vdw_geometric(x_aii, y_aii, x_bii, y_bii, r6, &v_a, &v_b);
    } else {
        calc_vdw_arithmetic(x_aii, y_aii, x_bii, y_bii, r6, &v_a, &v_b);
    }
    v_a *= lambda;
    v_b *= lambda;
    evdw = v_a - v_b;
    dv = r2 * (-ecoul - 12.0 * v_a + 6.0 * v_b);
}

__global__ void calc_nonbonded_force_kernel(
    const int nx,
    const int ny,

    const int* x_charges_types,
    const int* y_charges_types,
    const ccharge_t* ccharges_table,

    const int* x_atypes_types,
    const int* y_atypes_types,
    const catype_t* catypes_table,

    const topo_t d_topo,

    const bool* d_excluded,

    const int* d_LJ_matrix,

    const int* x_idx_list,
    const int* y_idx_list,

    const coord_t* d_coords,

    dvel_t* d_dvelocities,

    double* evdw_tot,
    double* ecoul_tot,

    bool symmetric,

    bool disable_water_h_lj,

    // helper variables
    const int n_atoms_solute,
    const int n_qelscales,
    const double lambda,
    const q_elscale_t* d_qelscales  // todo: Now doesn't use it. Should optimize it later

) {
    const int x_block_num = (nx + 31) >> 5;
    const int y_block_num = (ny + 31) >> 5;

    const int total_tiles = symmetric ? (x_block_num * (x_block_num + 1)) >> 1 : x_block_num * y_block_num;

    const int warps_per_block = blockDim.x >> 5;
    const int tid = threadIdx.x;
    const int lane = tid & 31;
    const int warp_in_block = tid >> 5;

    const int tile = blockIdx.x * warps_per_block + warp_in_block;
    if (tile >= total_tiles) return;

    int tile_x, tile_y;
    if (symmetric) {
        idx2xy(x_block_num, tile, tile_x, tile_y);
    } else {
        tile_x = tile % x_block_num;
        tile_y = tile / x_block_num;
    }

    const int base_x = tile_x << 5;
    const int base_y = tile_y << 5;

    int x_idx = base_x + lane;
    int y_idx = base_y + lane;

    int x_atom_idx = (x_idx < nx) ? x_idx_list[x_idx] : -1;
    int y_atom_idx = (y_idx < ny) ? y_idx_list[y_idx] : -1;

    coord_t invalid = {-1e9, -1e9, -1e9};
    coord_t x_coord = (x_atom_idx >= 0) ? d_coords[x_atom_idx] : invalid;
    coord_t y_coord = (y_atom_idx >= 0) ? d_coords[y_atom_idx] : invalid;

    bool x_excluded = (x_atom_idx >= 0) ? d_excluded[x_atom_idx] : true;
    bool y_excluded = (y_atom_idx >= 0) ? d_excluded[y_atom_idx] : true;

    int x_charge_type_idx = (x_idx < nx) ? x_charges_types[x_idx] : -1;
    int y_charge_type_idx = (y_idx < ny) ? y_charges_types[y_idx] : -1;
    double x_charge = (x_atom_idx >= 0 && x_charge_type_idx >= 0) ? ccharges_table[x_charge_type_idx].charge : 0.0;
    double y_charge = (y_atom_idx >= 0 && y_charge_type_idx >= 0) ? ccharges_table[y_charge_type_idx].charge : 0.0;

    int x_catype_type_idx = (x_idx < nx) ? x_atypes_types[x_idx] : -1;
    int y_catype_type_idx = (y_idx < ny) ? y_atypes_types[y_idx] : -1;
    catype_t x_type = (x_atom_idx >= 0 && x_catype_type_idx >= 0) ? catypes_table[x_catype_type_idx] : catype_t{};
    catype_t y_type = (y_atom_idx >= 0 && y_catype_type_idx >= 0) ? catypes_table[y_catype_type_idx] : catype_t{};

    double3 x_force = {0.0, 0.0, 0.0};
    double3 y_force = {0.0, 0.0, 0.0};

    double evdw_sum = 0.0;
    double ecoul_sum = 0.0;

    const unsigned mask = 0xffffffffu;

    auto lj_code = [&](int ix, int iy) -> int {
        if (ix < 0 || iy < 0 || ix >= n_atoms_solute || iy >= n_atoms_solute) return 0;
        return d_LJ_matrix[ix * n_atoms_solute + iy];
    };

    auto is_valid = [&]() -> bool {
        if (x_idx >= nx || y_idx >= ny) return false;

        if (symmetric && tile_x == tile_y && x_idx <= y_idx) return false;

        if (x_excluded || y_excluded) return false;

        // Skip interactions within the same rigid water molecule
        if (x_atom_idx >= n_atoms_solute && y_atom_idx >= n_atoms_solute) {
            int wx = (x_atom_idx - n_atoms_solute) / 3;
            int wy = (y_atom_idx - n_atoms_solute) / 3;
            if (wx == wy) return false;
        }

        bool bond23 = lj_code(x_atom_idx, y_atom_idx) == 3;
        if (bond23) return false;

        return true;
    };

    auto do_shuffle = [&]() {
        const int src = (lane + 1) & 31;
        y_idx = __shfl_sync(mask, y_idx, src);
        y_atom_idx = __shfl_sync(mask, y_atom_idx, src);
        y_coord = shfl_coord(y_coord, src, mask);
        y_excluded = __shfl_sync(mask, y_excluded, src);
        y_charge = __shfl_sync(mask, y_charge, src);
        y_type = shfl_catype(y_type, src, mask);

        y_force.x = shfl(y_force.x, src, mask);
        y_force.y = shfl(y_force.y, src, mask);
        y_force.z = shfl(y_force.z, src, mask);
    };

    if (disable_water_h_lj) {
        if (x_atom_idx >= n_atoms_solute && ((x_atom_idx - n_atoms_solute) % 3 != 0)) {
            x_type.aii_normal = 0.0;
            x_type.bii_normal = 0.0;
        }
        if (y_atom_idx >= n_atoms_solute && ((y_atom_idx - n_atoms_solute) % 3 != 0)) {
            y_type.aii_normal = 0.0;
            y_type.bii_normal = 0.0;
        }
    }

    for (int i = 0; i < 32; i++) {
        if (is_valid()) {
            bool bond14 = lj_code(x_atom_idx, y_atom_idx) == 1;
            double scaling = bond14 ? d_topo.el14_scale : 1.0;
            double ai_aii = bond14 ? x_type.aii_1_4 : x_type.aii_normal;
            double bi_bii = bond14 ? x_type.bii_1_4 : x_type.bii_normal;
            double aj_aii = bond14 ? y_type.aii_1_4 : y_type.aii_normal;
            double bj_bii = bond14 ? y_type.bii_1_4 : y_type.bii_normal;

            // todo: Now the idx is wrong, should optimize it later
            // for (int k = 0; k < n_qelscales; k++) {
            //     q_elscale_t qscale = d_qelscales[k];
            //     if ((x_charge_type_idx == qscale.qi) && (y_charge_type_idx == qscale.qj)) {
            //         scaling *= qscale.mu;
            //     }
            // }

            double evdw = 0, ecoul = 0, dv = 0;

            calculate_unforce_bound(
                x_coord,
                y_coord,
                x_charge,
                y_charge,
                ai_aii,
                aj_aii,
                bi_bii,
                bj_bii,
                d_topo.coulomb_constant,
                scaling,
                d_topo.vdw_rule,
                lambda,
                evdw,
                ecoul,
                dv);

            evdw_sum += evdw;
            ecoul_sum += ecoul;

            double3 d = {x_coord.x - y_coord.x, x_coord.y - y_coord.y, x_coord.z - y_coord.z};
            y_force.x -= dv * d.x;
            y_force.y -= dv * d.y;
            y_force.z -= dv * d.z;

            x_force.x += dv * d.x;
            x_force.y += dv * d.y;
            x_force.z += dv * d.z;
        }
        do_shuffle();
    }

    if (x_atom_idx >= 0) {
        atomicAdd(&d_dvelocities[x_atom_idx].x, x_force.x);
        atomicAdd(&d_dvelocities[x_atom_idx].y, x_force.y);
        atomicAdd(&d_dvelocities[x_atom_idx].z, x_force.z);
    }

    if (y_atom_idx >= 0) {
        atomicAdd(&d_dvelocities[y_atom_idx].x, y_force.x);
        atomicAdd(&d_dvelocities[y_atom_idx].y, y_force.y);
        atomicAdd(&d_dvelocities[y_atom_idx].z, y_force.z);
    }

    for (int offset = 16; offset > 0; offset >>= 1) {
        evdw_sum += __shfl_down_sync(mask, evdw_sum, offset);
        ecoul_sum += __shfl_down_sync(mask, ecoul_sum, offset);
    }
    if (lane == 0) {
        atomicAdd(evdw_tot, evdw_sum);
        atomicAdd(ecoul_tot, ecoul_sum);
    }
}

}  // namespace CudaNonbondedForce

std::pair<double, double> calc_nonbonded_force_host(
    int nx,
    int ny,
    int* x_idx_list,
    int* y_idx_list,
    bool symmetric,

    const int* x_charges_types,
    const int* y_charges_types,
    const ccharge_t* ccharges_table,

    const int* x_atypes_types,
    const int* y_atypes_types,
    const catype_t* catypes_table,
    const bool disable_water_h_lj, const double lambda) {
    using namespace CudaNonbondedForce;
    Context& host = Context::instance();
    CudaContext& context = CudaContext::instance();
    const int thread_num = 256;
    dim3 block_sz = dim3(thread_num);

    int tile_num_per_block = thread_num >> 5;

    int x_block_num = (nx + 31) >> 5;
    int y_block_num = (ny + 31) >> 5;

    int total_tiles = symmetric ? (x_block_num * (x_block_num + 1)) >> 1 : x_block_num * y_block_num;
    int grid_sz = (total_tiles + tile_num_per_block - 1) / tile_num_per_block;

    dim3 grid = dim3(grid_sz);

    cudaMemset(d_ecoul_total, 0, sizeof(double));
    cudaMemset(d_evdw_total, 0, sizeof(double));

    calc_nonbonded_force_kernel<<<grid, block_sz>>>(
        nx,
        ny,
        x_charges_types,
        y_charges_types,
        ccharges_table,
        x_atypes_types,
        y_atypes_types,
        catypes_table,
        host.topo,
        context.d_excluded,
        context.d_LJ_matrix,
        x_idx_list,
        y_idx_list,
        context.d_coords,
        context.d_dvelocities,
        d_evdw_total,
        d_ecoul_total,
        symmetric,
        disable_water_h_lj,
        host.n_atoms_solute,
        host.n_qelscales,
        lambda,
        context.d_q_elscales);

    cudaDeviceSynchronize();

    double evdw_tot = 0, ecoul_tot = 0;
    cudaMemcpy(&evdw_tot, d_evdw_total, sizeof(double), cudaMemcpyDeviceToHost);
    cudaMemcpy(&ecoul_tot, d_ecoul_total, sizeof(double), cudaMemcpyDeviceToHost);

    return {evdw_tot, ecoul_tot};
}

void init_nonbonded_force_kernel_data() {
    using namespace CudaNonbondedForce;
    if (!is_initialized) {
        check_cudaMalloc((void**)&d_evdw_total, sizeof(double));
        check_cudaMalloc((void**)&d_ecoul_total, sizeof(double));
        is_initialized = true;
    }
}

void cleanup_nonbonded_force() {
    using namespace CudaNonbondedForce;
    if (is_initialized) {
        cudaFree(d_evdw_total);
        cudaFree(d_ecoul_total);
        is_initialized = false;
    }
}
