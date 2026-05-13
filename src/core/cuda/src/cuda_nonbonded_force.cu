#include <iostream>

#include "context.h"
#include "cuda_nonbonded_force.cuh"
#include "constants.h"
#include "cuda_utility.cuh"

namespace CudaNonbondedForce {
bool is_initialized = false;
real_t *d_evdw_total, *d_ecoul_total;

template <typename WorkT>
struct nonbond_vec_t {
    WorkT x;
    WorkT y;
    WorkT z;
};

__device__ __forceinline__ float nonbond_rsqrt(float value) {
    return rsqrtf(value);
}

#ifndef QDYN_SPFP
__device__ __forceinline__ double nonbond_rsqrt(double value) {
    return rsqrt(value);
}
#endif

__device__ __forceinline__ void idx2xy(int n, int t, int& x, int& y) {
    x = (int)floorf((2 * n + 1 - sqrtf((2 * n + 1) * (2 * n + 1) - 8 * t)) * 0.5f);
    y = t - (x * n - (x * (x - 1) >> 1));
    if (y < 0) {
        x--;
        y += (n - x);
    }
    y += x;
}

template <typename T>
__device__ __forceinline__ T shfl_value(T v, int srcLane, unsigned mask = 0xffffffffu) {
    return __shfl_sync(mask, v, srcLane);
}

#ifndef QDYN_SPFP
template <>
__device__ __forceinline__ double shfl_value(double v, int srcLane, unsigned mask) {
    int2 a = *reinterpret_cast<int2*>(&v);
    a.x = __shfl_sync(mask, a.x, srcLane);
    a.y = __shfl_sync(mask, a.y, srcLane);
    return *reinterpret_cast<double*>(&a);
}
#endif

__device__ __forceinline__ coord_t shfl_coord(coord_t v, int srcLane, unsigned mask = 0xffffffffu) {
    v.x = shfl_value(v.x, srcLane, mask);
    v.y = shfl_value(v.y, srcLane, mask);
    v.z = shfl_value(v.z, srcLane, mask);
    return v;
}

template <typename WorkT>
__device__ void calculate_unforce_bound(
    const coord_t& x,
    const coord_t& y,

    const real_t charge_product,
    const vdw_pair_param_t& pair_param,

    const WorkT coulomb_constant,

    const WorkT scaling,
    const WorkT lambda,

    WorkT& evdw,
    WorkT& ecoul,
    WorkT& dv) {
    const WorkT dx = static_cast<WorkT>(x.x - y.x);
    const WorkT dy = static_cast<WorkT>(x.y - y.y);
    const WorkT dz = static_cast<WorkT>(x.z - y.z);
    const WorkT r = nonbond_rsqrt(dx * dx + dy * dy + dz * dz);
    const WorkT r2 = r * r;
    const WorkT r6 = r2 * r2 * r2;
    // real_t v_a = r6 * r6;
    // real_t v_b = r6;
    // ecoul = r;
    // evdw = v_a - v_b;
    // dv = r2 * (-ecoul - v_a + v_b);

    ecoul = scaling * coulomb_constant * charge_product * r * lambda;

    const WorkT v_a = static_cast<WorkT>(pair_param.a) * r6 * r6 * lambda;
    const WorkT v_b = static_cast<WorkT>(pair_param.b) * r6 * lambda;
    evdw = v_a - v_b;
    dv = r2 * (-ecoul - static_cast<WorkT>(12.0) * v_a + static_cast<WorkT>(6.0) * v_b);
}

template <typename WorkT>
__global__ void calc_nonbonded_force_kernel(
    const int nx,
    const int ny,

    const int* x_charges_types,
    const int* y_charges_types,
    const real_t* charge_pair_products,

    const int* x_atypes_types,
    const int* y_atypes_types,
    const vdw_pair_param_t* catype_pair_params,

    const topo_t d_topo,

    const bool* d_excluded,

    const int* d_LJ_matrix,

    const int* x_idx_list,
    const int* y_idx_list,

    const coord_t* d_coords,

    dvel_t* d_dvelocities,

    real_t* evdw_tot,
    real_t* ecoul_tot,

    bool symmetric,

    bool disable_water_h_lj,

    // helper variables
    const int n_atoms_solute,
    const int n_charge_types,
    const int zero_charge_type,
    const int n_catype_types,
    const int zero_catype_type,
    const int n_qelscales,
    const real_t lambda,
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

    coord_t invalid = {static_cast<real_t>(-1e9), static_cast<real_t>(-1e9), static_cast<real_t>(-1e9)};
    coord_t x_coord = (x_atom_idx >= 0) ? d_coords[x_atom_idx] : invalid;
    coord_t y_coord = (y_atom_idx >= 0) ? d_coords[y_atom_idx] : invalid;

    bool x_excluded = (x_atom_idx >= 0) ? d_excluded[x_atom_idx] : true;
    bool y_excluded = (y_atom_idx >= 0) ? d_excluded[y_atom_idx] : true;

    int x_charge_type_idx = (x_idx < nx) ? x_charges_types[x_idx] : zero_charge_type;
    int y_charge_type_idx = (y_idx < ny) ? y_charges_types[y_idx] : zero_charge_type;

    int x_catype_type_idx = (x_idx < nx) ? x_atypes_types[x_idx] : -1;
    int y_catype_type_idx = (y_idx < ny) ? y_atypes_types[y_idx] : -1;

    nonbond_vec_t<WorkT> x_force = {0.0, 0.0, 0.0};
    nonbond_vec_t<WorkT> y_force = {0.0, 0.0, 0.0};

    real_t evdw_sum = 0.0;
    real_t ecoul_sum = 0.0;

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

        int lj_relation = lj_code(x_atom_idx, y_atom_idx);
        if (lj_relation == 1 || lj_relation == 3) return false;

        return true;
    };

    auto do_shuffle = [&]() {
        const int src = (lane + 1) & 31;
        y_idx = __shfl_sync(mask, y_idx, src);
        y_atom_idx = __shfl_sync(mask, y_atom_idx, src);
        y_coord = shfl_coord(y_coord, src, mask);
        y_excluded = __shfl_sync(mask, y_excluded, src);
        y_charge_type_idx = __shfl_sync(mask, y_charge_type_idx, src);
        y_catype_type_idx = __shfl_sync(mask, y_catype_type_idx, src);

        y_force.x = shfl_value(y_force.x, src, mask);
        y_force.y = shfl_value(y_force.y, src, mask);
        y_force.z = shfl_value(y_force.z, src, mask);
    };

    if (disable_water_h_lj) {
        if (x_atom_idx >= n_atoms_solute && ((x_atom_idx - n_atoms_solute) % 3 != 0) && x_catype_type_idx >= 0) {
            x_catype_type_idx = zero_catype_type;
        }
        if (y_atom_idx >= n_atoms_solute && ((y_atom_idx - n_atoms_solute) % 3 != 0) && y_catype_type_idx >= 0) {
            y_catype_type_idx = zero_catype_type;
        }
    }

    const WorkT kernel_lambda = static_cast<WorkT>(lambda);
    const WorkT coulomb_constant = static_cast<WorkT>(d_topo.coulomb_constant);
    const int charge_pair_row = x_charge_type_idx * n_charge_types;
    const int pair_row = (x_catype_type_idx >= 0) ? x_catype_type_idx * n_catype_types : 0;

    for (int i = 0; i < 32; i++) {
        if (is_valid()) {
            WorkT scaling = static_cast<WorkT>(1.0);
            real_t charge_product = charge_pair_products[charge_pair_row + y_charge_type_idx];
            vdw_pair_param_t pair_param = catype_pair_params[pair_row + y_catype_type_idx];

            // todo: Now the idx is wrong, should optimize it later
            // for (int k = 0; k < n_qelscales; k++) {
            //     q_elscale_t qscale = d_qelscales[k];
            //     if ((x_charge_type_idx == qscale.qi) && (y_charge_type_idx == qscale.qj)) {
            //         scaling *= qscale.mu;
            //     }
            // }

            WorkT evdw = 0, ecoul = 0, dv = 0;

            calculate_unforce_bound(
                x_coord,
                y_coord,
                charge_product,
                pair_param,
                coulomb_constant,
                scaling,
                kernel_lambda,
                evdw,
                ecoul,
                dv);

            evdw_sum += evdw;
            ecoul_sum += ecoul;

            const WorkT dx = static_cast<WorkT>(x_coord.x - y_coord.x);
            const WorkT dy = static_cast<WorkT>(x_coord.y - y_coord.y);
            const WorkT dz = static_cast<WorkT>(x_coord.z - y_coord.z);
            y_force.x -= dv * dx;
            y_force.y -= dv * dy;
            y_force.z -= dv * dz;

            x_force.x += dv * dx;
            x_force.y += dv * dy;
            x_force.z += dv * dz;
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

std::pair<real_t, real_t> calc_nonbonded_force_host(
    int nx,
    int ny,
    int* x_idx_list,
    int* y_idx_list,
    bool symmetric,

    const int* x_charges_types,
    const int* y_charges_types,
    const int* x_atypes_types,
    const int* y_atypes_types,
    const bool disable_water_h_lj, const real_t lambda) {
    using namespace CudaNonbondedForce;
    Context& host = Context::instance();
    const int thread_num = 256;
    dim3 block_sz = dim3(thread_num);

    int tile_num_per_block = thread_num >> 5;

    int x_block_num = (nx + 31) >> 5;
    int y_block_num = (ny + 31) >> 5;

    int total_tiles = symmetric ? (x_block_num * (x_block_num + 1)) >> 1 : x_block_num * y_block_num;
    int grid_sz = (total_tiles + tile_num_per_block - 1) / tile_num_per_block;

    dim3 grid = dim3(grid_sz);

    cudaMemset(d_ecoul_total, 0, sizeof(real_t));
    cudaMemset(d_evdw_total, 0, sizeof(real_t));

    auto launch_kernel = [&](auto work_tag) {
        using WorkT = decltype(work_tag);
        calc_nonbonded_force_kernel<WorkT><<<grid, block_sz>>>(
            nx,
            ny,
            x_charges_types,
            y_charges_types,
            host.charge_pair_products->gpu_data_p,
            x_atypes_types,
            y_atypes_types,
            host.catype_pair_params->gpu_data_p,
            host.topo,
            host.excluded->gpu_data_p,
            host.LJ_matrix->gpu_data_p,
            x_idx_list,
            y_idx_list,
            host.coords->gpu_data_p,
            host.dvelocities->gpu_data_p,
            d_evdw_total,
            d_ecoul_total,
            symmetric,
            disable_water_h_lj,
            host.n_atoms_solute,
            host.n_charge_types,
            host.zero_charge_type,
            host.n_catype_types,
            host.zero_catype_type,
            host.n_qelscales,
            lambda,
            host.q_elscales->gpu_data_p);
    };

    launch_kernel(nonbond_work_t{});

    cudaDeviceSynchronize();

    real_t evdw_tot = 0, ecoul_tot = 0;
    cudaMemcpy(&evdw_tot, d_evdw_total, sizeof(real_t), cudaMemcpyDeviceToHost);
    cudaMemcpy(&ecoul_tot, d_ecoul_total, sizeof(real_t), cudaMemcpyDeviceToHost);

    return {evdw_tot, ecoul_tot};
}

void init_nonbonded_force_kernel_data() {
    using namespace CudaNonbondedForce;
    if (!is_initialized) {
        check_cudaMalloc((void**)&d_evdw_total, sizeof(real_t));
        check_cudaMalloc((void**)&d_ecoul_total, sizeof(real_t));
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
