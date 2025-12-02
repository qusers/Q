#pragma once
#include "cuda/include/CudaContext.cuh"
#include "cuda/include/CudaNonbondedWWForce.cuh"
namespace CudaNonbondedWWForce {

bool is_initialized = false;
double *D_WW_evdw_TOT, *D_WW_ecoul_TOT, WW_evdw_TOT, WW_ecoul_TOT;

__device__ __forceinline__ void calculate_unforce_bound(
    const int y, const int x, const coord_t& q, const coord_t& p,
    const topo_t& D_topo, const double& crg_ow, const double& crg_hw,
    const double& A_OO, const double& B_OO, double& evdw, double& ecoul,
    double& dv, double& tmpx, double& tmpy, double& tmpz) {
    int belong_y = y / 3;
    int belong_x = x / 3;
    if (belong_y == belong_x) {
        return;
    }

    bool y_is_o = (y % 3 == 0);
    bool x_is_o = (x % 3 == 0);

    // Compute distance components
    tmpx = p.x - q.x;
    tmpy = p.y - q.y;
    tmpz = p.z - q.z;
    // double inv_dis = 1.0 / sqrt(pow(tmpx, 2) + pow(tmpy, 2) + pow(tmpz, 2));
    double inv_dis = rsqrt(tmpx * tmpx + tmpy * tmpy + tmpz * tmpz);
    double inv_dis2 = inv_dis * inv_dis;

    ecoul = inv_dis * D_topo.coulomb_constant * (y_is_o ? crg_ow : crg_hw) *
            (x_is_o ? crg_ow : crg_hw);
    double v_a = 0, v_b = 0;
    if (y_is_o && x_is_o) {
        double inv_dis6 = inv_dis2 * inv_dis2 * inv_dis2;
        double inv_dis12 = inv_dis6 * inv_dis6;
        v_a = A_OO * inv_dis12;
        v_b = B_OO * inv_dis6;
        evdw = v_a - v_b;
        dv = inv_dis * inv_dis * (-ecoul - 12.0 * v_a + 6.0 * v_b);
    } else {
        dv = inv_dis * inv_dis * -ecoul;
    }
}

template <const int Thread_x, const int Thread_y, const int Block_x,
          const int Block_y>
__global__ void
calc_ww(const int N, const double crg_ow, const double crg_hw,
        const double A_OO, const double B_OO, const topo_t D_topo,
        coord_t* __restrict__ W, dvel_t* __restrict__ DV_W,
        double* __restrict__ Evdw_TOT, double* __restrict__ ecoul_TOT) {
    // Calculate block boundaries
    int NX = N;
    int NY = (N + 1) / 2;
    int x_cal_num = blockDim.x * Block_x;
    int y_cal_num = blockDim.y * Block_y;
    int block_x_left_begin = 1 + blockIdx.x * x_cal_num;
    int block_y_left_begin = blockIdx.y * y_cal_num;
    int block_x_left_end = min(block_x_left_begin + x_cal_num - 1, NX - 1);
    int block_y_left_end = min(block_y_left_begin + y_cal_num - 1, NY - 1);

    // Shared memory declarations with padding to avoid bank conflicts
    __shared__ coord_t p[2 * Thread_x * Block_x + 1];
    __shared__ coord_t q[2 * Thread_y * Block_y + 1];
    __shared__ double sum_row_x[2 * Thread_y * Block_y + 1];
    __shared__ double sum_row_y[2 * Thread_y * Block_y + 1];
    __shared__ double sum_row_z[2 * Thread_y * Block_y + 1];

    __shared__ double block_ecoul[(Thread_x * Thread_y + 31) / 32],
        block_evdw[(Thread_x * Thread_y + 31) / 32];

    // Thread indices
    int thread_y_left_begin = block_y_left_begin + threadIdx.y * Block_y;
    int thread_x_left_begin = block_x_left_begin + threadIdx.x * Block_x;
    int cur_thread_num = blockDim.x * threadIdx.y + threadIdx.x;
    int thread_num = blockDim.x * blockDim.y;

// Optimized coordinate loading with coalesced memory access
#pragma unroll
    for (int i = cur_thread_num; i < x_cal_num && block_x_left_begin + i < NX; i += thread_num) {
        p[i] = W[block_x_left_begin + i];
        p[x_cal_num + i] = W[N - (block_x_left_begin + i)];
    }

#pragma unroll
    for (int i = cur_thread_num; i < y_cal_num && block_y_left_begin + i < NY; i += thread_num) {
        q[i] = W[block_y_left_begin + i];
        q[y_cal_num + i] = W[N - 1 - (block_y_left_begin + i)];
    }

// Initialize sum arrays
#pragma unroll
    for (int i = cur_thread_num; i < y_cal_num && block_y_left_begin + i < NY; i += thread_num) {
        sum_row_x[i] = sum_row_x[y_cal_num + i] = 0;
        sum_row_y[i] = sum_row_y[y_cal_num + i] = 0;
        sum_row_z[i] = sum_row_z[y_cal_num + i] = 0;
    }
    __syncthreads();
    // Initialize column sums
    double sum_col_x[2 * Block_x], sum_col_y[2 * Block_x], sum_col_z[2 * Block_x];
#pragma unroll
    for (int i = 0; i < Block_x; i++) {
        sum_col_x[i] = sum_col_x[Block_x + i] = 0;
        sum_col_y[i] = sum_col_y[Block_x + i] = 0;
        sum_col_z[i] = sum_col_z[Block_x + i] = 0;
    }

    // Main computation loop with reduced thread divergence
    double evdw_tot = 0, ecoul_tot = 0;
#pragma unroll
    for (int i = 0; i < Block_x; i++) {
        int i2 = i;
        int x = thread_x_left_begin + i2;
        if (x >= NX) continue;

        int offset_x = x - block_x_left_begin;
        const coord_t& now_p0 = p[offset_x];
        const coord_t& now_p1 = p[x_cal_num + offset_x];

#pragma unroll
        for (int j = 0; j < Block_y; j++) {
            int j2 = (j + threadIdx.x) % Block_y;
            int y = thread_y_left_begin + j2;
            if (y >= NY) continue;

            // Optimized condition check
            if (y >= x && (N % 2 == 0 || y != NY - 1)) {
                int offset_y = y - block_y_left_begin;
                double evdw = 0, ecoul = 0, dv = 0;
                double tmpx = 0, tmpy = 0, tmpz = 0;

                int y2 = N - 1 - y;
                int x2 = N - x;
                calculate_unforce_bound(y2, x2, q[y_cal_num + offset_y], now_p1, D_topo,
                                        crg_ow, crg_hw, A_OO, B_OO, evdw, ecoul, dv,
                                        tmpx, tmpy, tmpz);

                evdw_tot += evdw;
                ecoul_tot += ecoul;
                double v_x = dv * tmpx;
                double v_y = dv * tmpy;
                double v_z = dv * tmpz;

                sum_col_x[Block_x + i2] += v_x;
                sum_col_y[Block_x + i2] += v_y;
                sum_col_z[Block_x + i2] += v_z;
                sum_row_x[y_cal_num + offset_y] -= v_x;
                sum_row_y[y_cal_num + offset_y] -= v_y;
                sum_row_z[y_cal_num + offset_y] -= v_z;
            } else if (y < x) {
                int offset_y = y - block_y_left_begin;
                double evdw = 0, ecoul = 0, dv = 0;
                double tmpx = 0, tmpy = 0, tmpz = 0;

                calculate_unforce_bound(y, x, q[offset_y], now_p0, D_topo, crg_ow,
                                        crg_hw, A_OO, B_OO, evdw, ecoul, dv, tmpx, tmpy,
                                        tmpz);

                evdw_tot += evdw;
                ecoul_tot += ecoul;
                double v_x = dv * tmpx;
                double v_y = dv * tmpy;
                double v_z = dv * tmpz;

                sum_col_x[i2] += v_x;
                sum_col_y[i2] += v_y;
                sum_col_z[i2] += v_z;
                sum_row_x[offset_y] -= v_x;
                sum_row_y[offset_y] -= v_y;
                sum_row_z[offset_y] -= v_z;
            }
        }
    }

// Optimized reduction using warp-level primitives
#pragma unroll
    for (unsigned w = 16; w >= 1; w /= 2) {
        ecoul_tot += __shfl_down_sync(0xffffffff, ecoul_tot, w);
        evdw_tot += __shfl_down_sync(0xffffffff, evdw_tot, w);
    }

    // Store block results
    int warp_id = cur_thread_num / warpSize;
    int lane_id = cur_thread_num % warpSize;
    if (lane_id == 0) {
        block_evdw[warp_id] = evdw_tot;
        block_ecoul[warp_id] = ecoul_tot;
    }
    __syncthreads();

    // Final reduction
    if (warp_id == 0) {
        double val1 = (lane_id < (Thread_x * Thread_y + 31) / 32) ? block_ecoul[lane_id] : 0;
        double val2 = (lane_id < (Thread_x * Thread_y + 31) / 32) ? block_evdw[lane_id] : 0;

#pragma unroll
        for (unsigned w = 16; w >= 1; w /= 2) {
            val1 += __shfl_down_sync(0xffffffff, val1, w);
            val2 += __shfl_down_sync(0xffffffff, val2, w);
        }

        if (lane_id == 0) {
            atomicAdd(ecoul_TOT, val1);
            atomicAdd(Evdw_TOT, val2);
        }
    }

// Optimized row reduction
#pragma unroll
    for (int i = cur_thread_num; i < y_cal_num && block_y_left_begin + i < NY; i += thread_num) {
        int idx = block_y_left_begin + i;
        if (idx < block_x_left_end) {
            atomicAdd(&DV_W[idx].x, sum_row_x[i]);
            atomicAdd(&DV_W[idx].y, sum_row_y[i]);
            atomicAdd(&DV_W[idx].z, sum_row_z[i]);
        }
        if (idx >= block_x_left_begin) {
            idx = N - 1 - idx;
            atomicAdd(&DV_W[idx].x, sum_row_x[i + y_cal_num]);
            atomicAdd(&DV_W[idx].y, sum_row_y[i + y_cal_num]);
            atomicAdd(&DV_W[idx].z, sum_row_z[i + y_cal_num]);
        }
    }

// Optimized column reduction
#pragma unroll
    for (int i = 0; i < Block_x; i++) {
        int i2 = (i + threadIdx.y) % Block_x;
        int idx = thread_x_left_begin + i2;
        if (idx < N) {
            if (block_y_left_begin < idx) {
                atomicAdd(&DV_W[idx].x, sum_col_x[i2]);
                atomicAdd(&DV_W[idx].y, sum_col_y[i2]);
                atomicAdd(&DV_W[idx].z, sum_col_z[i2]);
            }
            if (block_y_left_end >= idx) {
                atomicAdd(&DV_W[N - idx].x, sum_col_x[i2 + Block_x]);
                atomicAdd(&DV_W[N - idx].y, sum_col_y[i2 + Block_x]);
                atomicAdd(&DV_W[N - idx].z, sum_col_z[i2 + Block_x]);
            }
        }
    }
}
}  // namespace CudaNonbondedWWForce
void calc_nonbonded_ww_forces_host_v2() {
    using namespace CudaNonbondedWWForce;
    int N = 3 * n_waters;
    int mem_size_W = N * sizeof(coord_t);
    int mem_size_DV_W = N * sizeof(dvel_t);

    WW_evdw_TOT = 0;
    WW_ecoul_TOT = 0;
    cudaMemcpy(D_WW_evdw_TOT, &WW_evdw_TOT, sizeof(double),
               cudaMemcpyHostToDevice);
    cudaMemcpy(D_WW_ecoul_TOT, &WW_ecoul_TOT, sizeof(double),
               cudaMemcpyHostToDevice);

    CudaContext& ctx = CudaContext::instance();
    auto W = ctx.d_coords + n_atoms_solute;
    auto DV_W = ctx.d_dvelocities + n_atoms_solute;

    // Optimize thread block configuration
    const int thread_num_x = 32;  // Keep at 32 for better warp utilization
    const int thread_num_y = 1;   // Keep at 1 for better memory coalescing
    dim3 block_sz = dim3(thread_num_x, thread_num_y);

    const int N_ITERATION_Y = 32;  // Keep at 32 for better memory access pattern
    const int N_ITERATION_X = 1;   // Keep at 1 for better memory coalescing

    int block_element_sz_x = N_ITERATION_X * block_sz.x;
    int block_element_sz_y = N_ITERATION_Y * block_sz.y;

    int grid_sz_x = (N - 1 + block_element_sz_x - 1) / (block_element_sz_x);
    int grid_sz_y = ((N + 1) / 2 + block_element_sz_y - 1) / block_element_sz_y;
    dim3 grid = dim3(grid_sz_x, grid_sz_y);

    calc_ww<thread_num_x, thread_num_y, N_ITERATION_X, N_ITERATION_Y>
        <<<grid, block_sz>>>(N, crg_ow, crg_hw, A_OO, B_OO, topo, W, DV_W,
                             D_WW_evdw_TOT, D_WW_ecoul_TOT);

    cudaDeviceSynchronize();
    cudaMemcpy(&dvelocities[n_atoms_solute], DV_W, mem_size_DV_W,
               cudaMemcpyDeviceToHost);
    cudaMemcpy(&WW_evdw_TOT, D_WW_evdw_TOT, sizeof(double),
               cudaMemcpyDeviceToHost);
    cudaMemcpy(&WW_ecoul_TOT, D_WW_ecoul_TOT, sizeof(double),
               cudaMemcpyDeviceToHost);
    E_nonbond_ww.Uvdw += WW_evdw_TOT;
    E_nonbond_ww.Ucoul += WW_ecoul_TOT;
}

void init_nonbonded_ww_force_kernel_data() {
    using namespace CudaNonbondedWWForce;
    if (!is_initialized) {
        catype_t catype_ow;                // Atom type of first O, H atom
        ccharge_t ccharge_ow, ccharge_hw;  // Charge of first O, H atom

        catype_ow = catypes[atypes[n_atoms_solute].code - 1];
        ccharge_ow = ccharges[charges[n_atoms_solute].code - 1];
        ccharge_hw = ccharges[charges[n_atoms_solute + 1].code - 1];

        A_OO = pow(catype_ow.aii_normal, 2);
        B_OO = pow(catype_ow.bii_normal, 2);

        crg_ow = ccharge_ow.charge;
        crg_hw = ccharge_hw.charge;

        check_cudaMalloc((void**)&D_WW_evdw_TOT, sizeof(double));
        check_cudaMalloc((void**)&D_WW_ecoul_TOT, sizeof(double));
        is_initialized = true;
    }
}

void cleanup_nonbonded_ww_force() {
    using namespace CudaNonbondedWWForce;
    if (is_initialized) {
        cudaFree(D_WW_evdw_TOT);
        cudaFree(D_WW_ecoul_TOT);
        is_initialized = false;
    }
}