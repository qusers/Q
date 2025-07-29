#include "solvent.h"
#include "system.h"
#include "utils.h"

#include <stdio.h>
#include <time.h>
#include <unistd.h>

/* =============================================
 * == SOLVENT INTERACTIONS
 * =============================================
 */

coord_t *W;
dvel_t *DV_W;
calc_ww_t *WW_MAT, *h_WW_MAT;
calc_pw_t *PW_MAT, *h_PW_MAT;

double *D_WW_Evdw, *D_WW_Ecoul, *h_WW_Evdw, *h_WW_Ecoul;
double *D_WW_evdw_TOT, *D_WW_ecoul_TOT, WW_evdw_TOT, WW_ecoul_TOT;

double *D_PW_Evdw, *D_PW_Ecoul, *h_PW_Evdw, *h_PW_Ecoul;
double *D_PW_evdw_TOT, *D_PW_ecoul_TOT, PW_evdw_TOT, PW_ecoul_TOT;

// Constants pointers
bool pw_gpu_set = false;

// ONLY call if there are actually solvent atoms, or get segfaulted
void calc_nonbonded_ww_forces() {
    double rOX, rH1X, rH2X, r2;
    coord_t dOX, dH1X, dH2X;
    double Vel, V_a, V_b, dv;

    // Initialize water constants
    if (A_OO == 0) {
        catype_t catype_ow; // Atom type of first O, H atom
        ccharge_t ccharge_ow, ccharge_hw; // Charge of first O, H atom

        catype_ow = catypes[atypes[n_atoms_solute].code - 1];
        ccharge_ow = ccharges[charges[n_atoms_solute].code - 1];
        ccharge_hw = ccharges[charges[n_atoms_solute + 1].code - 1];

        A_OO = pow(catype_ow.aii_normal, 2);
        B_OO = pow(catype_ow.bii_normal, 2);

        crg_ow = ccharge_ow.charge;
        crg_hw = ccharge_hw.charge;
    }

    // double ecoul_block = 0, evdw_block = 0;

    for (int i = n_atoms_solute; i < n_atoms; i += 3) {
        for (int j = i + 3; j < n_atoms; j += 3) {
            double ecoul = 0, evdw = 0;
            // if (excluded[i] || excluded[j]) continue;
            // --- O - (O,H1,H2) ---
            dOX.x = coords[j].x - coords[i].x;
            dOX.y = coords[j].y - coords[i].y;
            dOX.z = coords[j].z - coords[i].z;
            rOX = pow(dOX.x, 2) + pow(dOX.y, 2) + pow(dOX.z, 2);

            // O - H1
            j += 1;

            dH1X.x = coords[j].x - coords[i].x;
            dH1X.y = coords[j].y - coords[i].y;
            dH1X.z = coords[j].z - coords[i].z;
            rH1X = pow(dH1X.x, 2) + pow(dH1X.y, 2) + pow(dH1X.z, 2);

            // O-H2 (X=H2)
            j += 1;

            dH2X.x = coords[j].x - coords[i].x;
            dH2X.y = coords[j].y - coords[i].y;
            dH2X.z = coords[j].z - coords[i].z;
            rH2X = pow(dH2X.x, 2) + pow(dH2X.y, 2) + pow(dH2X.z, 2);
            rOX = sqrt(1 / rOX);
            rH1X = sqrt(1 / rH1X);
            rH2X = sqrt(1 / rH2X);

            // O - O
            r2 = rOX * rOX;
            Vel = topo.coulomb_constant * pow(crg_ow, 2) * rOX;
            V_a = A_OO * (r2 * r2 * r2) * (r2 * r2 * r2);
            V_b = B_OO * (r2 * r2 * r2);
            evdw += (V_a - V_b);
            ecoul += Vel;
            dv = r2 * (-Vel - 12 * V_a + 6 * V_b);
            j -= 2; // move pointer back to O in interacting molecule
            dvelocities[i].x -= (dv * dOX.x);
            dvelocities[j].x += (dv * dOX.x);
            dvelocities[i].y -= (dv * dOX.y);
            dvelocities[j].y += (dv * dOX.y);
            dvelocities[i].z -= (dv * dOX.z);
            dvelocities[j].z += (dv * dOX.z);

            // O - H1
            r2 = pow(rH1X, 2);
            Vel = topo.coulomb_constant * crg_ow * crg_hw * rH1X;
            ecoul += Vel;
            dv = r2 * (-Vel);
            j += 1; // point to H1 in j-molecule

            dvelocities[i].x -= (dv * dH1X.x);
            dvelocities[j].x += (dv * dH1X.x);
            dvelocities[i].y -= (dv * dH1X.y);
            dvelocities[j].y += (dv * dH1X.y);
            dvelocities[i].z -= (dv * dH1X.z);
            dvelocities[j].z += (dv * dH1X.z);

            // O - H2
            r2 = pow(rH2X, 2);
            Vel = topo.coulomb_constant * crg_ow * crg_hw * rH2X;
            ecoul += Vel;
            dv = r2 * (-Vel);
            j += 1; // point to H2 in j-molecule

            dvelocities[i].x -= (dv * dH2X.x);
            dvelocities[j].x += (dv * dH2X.x);
            dvelocities[i].y -= (dv * dH2X.y);
            dvelocities[j].y += (dv * dH2X.y);
            dvelocities[i].z -= (dv * dH2X.z);
            dvelocities[j].z += (dv * dH2X.z);

            // --- H1 - (O,H1,H2) ---
            i += 1; // Point to H1 in i-molecule
            j -= 2; // Point to O in j-molecule

            // H1 - O (X=O)
            dOX.x = coords[j].x - coords[i].x;
            dOX.y = coords[j].y - coords[i].y;
            dOX.z = coords[j].z - coords[i].z;
            rOX = pow(dOX.x, 2) + pow(dOX.y, 2) + pow(dOX.z, 2);

            // H1 - H1 (X=H1)
            j += 1;

            dH1X.x = coords[j].x - coords[i].x;
            dH1X.y = coords[j].y - coords[i].y;
            dH1X.z = coords[j].z - coords[i].z;
            rH1X = pow(dH1X.x, 2) + pow(dH1X.y, 2) + pow(dH1X.z, 2);

            // H1 - H2 (X=H2)
            j += 1;

            dH2X.x = coords[j].x - coords[i].x;
            dH2X.y = coords[j].y - coords[i].y;
            dH2X.z = coords[j].z - coords[i].z;
            rH2X = pow(dH2X.x, 2) + pow(dH2X.y, 2) + pow(dH2X.z, 2);
            rOX = sqrt(1 / rOX);
            rH1X = sqrt(1 / rH1X);
            rH2X = sqrt(1 / rH2X);

            // H1 - O
            r2 = rOX * rOX;
            Vel = topo.coulomb_constant * crg_hw * crg_ow * rOX;
            ecoul += Vel;
            dv = r2 * (-Vel);
            j -= 2; // move pointer back to O in interacting molecule
            dvelocities[i].x -= (dv * dOX.x);
            dvelocities[j].x += (dv * dOX.x);
            dvelocities[i].y -= (dv * dOX.y);
            dvelocities[j].y += (dv * dOX.y);
            dvelocities[i].z -= (dv * dOX.z);
            dvelocities[j].z += (dv * dOX.z);

            // H1 - H1
            r2 = pow(rH1X, 2);
            Vel = topo.coulomb_constant * crg_hw * crg_hw * rH1X;
            ecoul += Vel;
            dv = r2 * (-Vel);
            j += 1; // point to H1 in j-molecule
            dvelocities[i].x -= (dv * dH1X.x);
            dvelocities[j].x += (dv * dH1X.x);
            dvelocities[i].y -= (dv * dH1X.y);
            dvelocities[j].y += (dv * dH1X.y);
            dvelocities[i].z -= (dv * dH1X.z);
            dvelocities[j].z += (dv * dH1X.z);

            // H1 - H2
            r2 = pow(rH2X, 2);
            Vel = topo.coulomb_constant * crg_hw * crg_hw * rH2X;
            ecoul += Vel;
            dv = r2 * (-Vel);
            j += 1; // point to H2 in j-molecule
            dvelocities[i].x -= (dv * dH2X.x);
            dvelocities[j].x += (dv * dH2X.x);
            dvelocities[i].y -= (dv * dH2X.y);
            dvelocities[j].y += (dv * dH2X.y);
            dvelocities[i].z -= (dv * dH2X.z);
            dvelocities[j].z += (dv * dH2X.z);

            // --- H2 - (O,H1,H2) ---
            i += 1; // Point to H2 in i-molecule
            j -= 2; // Point to O in j-molecule

            // H2 - O (X=O)
            dOX.x = coords[j].x - coords[i].x;
            dOX.y = coords[j].y - coords[i].y;
            dOX.z = coords[j].z - coords[i].z;
            rOX = pow(dOX.x, 2) + pow(dOX.y, 2) + pow(dOX.z, 2);

            // H2 - H1 (X=H1)
            j += 1;

            dH1X.x = coords[j].x - coords[i].x;
            dH1X.y = coords[j].y - coords[i].y;
            dH1X.z = coords[j].z - coords[i].z;
            rH1X = pow(dH1X.x, 2) + pow(dH1X.y, 2) + pow(dH1X.z, 2);

            // H2 - H2 (X=H2)
            j += 1;

            dH2X.x = coords[j].x - coords[i].x;
            dH2X.y = coords[j].y - coords[i].y;
            dH2X.z = coords[j].z - coords[i].z;
            rH2X = pow(dH2X.x, 2) + pow(dH2X.y, 2) + pow(dH2X.z, 2);
            rOX = sqrt(1 / rOX);
            rH1X = sqrt(1 / rH1X);
            rH2X = sqrt(1 / rH2X);

            // H2 - O
            r2 = rOX * rOX;
            Vel = topo.coulomb_constant * crg_hw * crg_ow * rOX;
            ecoul += Vel;
            dv = r2 * (-Vel);
            j -= 2; // move pointer back to O in interacting molecule
            dvelocities[i].x -= (dv * dOX.x);
            dvelocities[j].x += (dv * dOX.x);
            dvelocities[i].y -= (dv * dOX.y);
            dvelocities[j].y += (dv * dOX.y);
            dvelocities[i].z -= (dv * dOX.z);
            dvelocities[j].z += (dv * dOX.z);

            // H2 - H1
            r2 = pow(rH1X, 2);
            Vel = topo.coulomb_constant * crg_hw * crg_hw * rH1X;
            ecoul += Vel;
            dv = r2 * (-Vel);
            j += 1; // point to H1 in j-molecule
            dvelocities[i].x -= (dv * dH1X.x);
            dvelocities[j].x += (dv * dH1X.x);
            dvelocities[i].y -= (dv * dH1X.y);
            dvelocities[j].y += (dv * dH1X.y);
            dvelocities[i].z -= (dv * dH1X.z);
            dvelocities[j].z += (dv * dH1X.z);

            // H1 - H2
            r2 = pow(rH2X, 2);
            Vel = topo.coulomb_constant * crg_hw * crg_hw * rH2X;
            ecoul += Vel;
            dv = r2 * (-Vel);
            j += 1; // point to H2 in j-molecule
            dvelocities[i].x -= (dv * dH2X.x);
            dvelocities[j].x += (dv * dH2X.x);
            dvelocities[i].y -= (dv * dH2X.y);
            dvelocities[j].y += (dv * dH2X.y);
            dvelocities[i].z -= (dv * dH2X.z);
            dvelocities[j].z += (dv * dH2X.z);

            // Put i and j indexes back to prevent troubles
            i -= 2;
            j -= 2;

            E_nonbond_ww.Uvdw += evdw;
            E_nonbond_ww.Ucoul += ecoul;
            // if (70*8 <= i/3 && i/3 < 71*8 && 70*8 <= j/3 && j/3 < 71*8) {
            //     ecoul_block += ecoul;
            //     evdw_block += evdw;
            // printf("Evdw[%d][%d] = %f\n", i/3, j/3, evdw);
            // printf("Ecoul[%d][%d] = %f\n", i/3, j/3, ecoul);
            // }
        }
    }

    // for (int i = 0; i < n_atoms; i++) {
    //     printf("dvelocities[%d] = %f %f %f\n", i, dvelocities[i].x,
    //     dvelocities[i].y, dvelocities[i].z);
    // }

    // printf("ecoul_block = %f, evdw_block = %f\n", ecoul_block, evdw_block);

#ifdef DEBUG
  printf("solvent: Ecoul = %f Evdw = %f\n", energies.Ucoul, energies.Uvdw);
#endif
}

//__device__ coord_t get_coord(int y, int x, int flag, int N,
//                             coord_t *__restrict__ W) {
//  if (flag == 0) {
//    if (y >= x) {
//      return W[N - x];
//    } else {
//      return W[x];
//    }
//  } else {
//    if (y >= x) {
//      return W[N - 1 - y];
//    } else {
//      return W[y];
//    }
//  }
//}

__device__ __forceinline__ void calculate_unforce_bound(
    const int y, const int x, const coord_t &q, const coord_t &p,
    const topo_t &D_topo, const double &crg_ow, const double &crg_hw,
    const double &A_OO, const double &B_OO, double &evdw, double &ecoul,
    double &dv, double &tmpx, double &tmpy, double &tmpz) {
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

template<const int Thread_x, const int Thread_y, const int Block_x,
    const int Block_y>
__global__ void
calc_ww(const int N, const double crg_ow, const double crg_hw,
        const double A_OO, const double B_OO, const topo_t D_topo,
        coord_t *__restrict__ W, dvel_t *__restrict__ DV_W,
        double *__restrict__ Evdw_TOT, double *__restrict__ ecoul_TOT) {
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
        const coord_t &now_p0 = p[offset_x];
        const coord_t &now_p1 = p[x_cal_num + offset_x];

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

void calc_nonbonded_ww_forces_host() {
    /*
     * mem_size_w: The number of water atom = 3 * n_waters
    W: The coords of the water, size: mem_size_w
    WW_MAT: The distance between each water atom
    DV_W: dvel of each water atom
    D_WW_Evdw:
    D_WW_Ecoul:

    D_WW_evdw_TOT:
    D_WW_ecoul_TOT:

    h_WW_MAT:
    h_WW_Evdw:
    h_WW_Ecoul:


    */

    int N = 3 * n_waters;
    int mem_size_W = N * sizeof(coord_t);
    int mem_size_DV_W = N * sizeof(dvel_t);

    if (A_OO == 0) {
        catype_t catype_ow; // Atom type of first O, H atom
        ccharge_t ccharge_ow, ccharge_hw; // Charge of first O, H atom

        catype_ow = catypes[atypes[n_atoms_solute].code - 1];
        ccharge_ow = ccharges[charges[n_atoms_solute].code - 1];
        ccharge_hw = ccharges[charges[n_atoms_solute + 1].code - 1];

        A_OO = pow(catype_ow.aii_normal, 2);
        B_OO = pow(catype_ow.bii_normal, 2);

        crg_ow = ccharge_ow.charge;
        crg_hw = ccharge_hw.charge;

        check_cudaMalloc((void **) &W, mem_size_W);
        check_cudaMalloc((void **) &DV_W, mem_size_DV_W);
        check_cudaMalloc((void **) &D_WW_evdw_TOT, sizeof(double));
        check_cudaMalloc((void **) &D_WW_ecoul_TOT, sizeof(double));
    }

    WW_evdw_TOT = 0;
    WW_ecoul_TOT = 0;
    cudaMemcpy(W, &coords[n_atoms_solute], mem_size_W, cudaMemcpyHostToDevice);
    cudaMemcpy(DV_W, &dvelocities[n_atoms_solute], mem_size_DV_W,
               cudaMemcpyHostToDevice);
    cudaMemcpy(D_WW_evdw_TOT, &WW_evdw_TOT, sizeof(double),
               cudaMemcpyHostToDevice);
    cudaMemcpy(D_WW_ecoul_TOT, &WW_ecoul_TOT, sizeof(double),
               cudaMemcpyHostToDevice);

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

void calc_nonbonded_pw_forces() {
    coord_t da;
    double r2a, ra, r6a;
    double Vela, V_a, V_b;
    double dva;
    double qi, qj;
    double ai_aii, aj_aii, ai_bii, aj_bii;
    catype_t ai_type, aj_type;
    int i;

    for (int pi = 0; pi < n_patoms; pi++) {
        for (int j = n_atoms_solute; j < n_atoms; j++) {
            i = p_atoms[pi].a - 1;
            if (excluded[i] || excluded[j])
                continue;
            qi = ccharges[charges[i].code - 1].charge;
            qj = ccharges[charges[j].code - 1].charge;

            ai_type = catypes[atypes[i].code - 1];
            aj_type = catypes[atypes[j].code - 1];

            da.x = coords[j].x - coords[i].x;
            da.y = coords[j].y - coords[i].y;
            da.z = coords[j].z - coords[i].z;
            r2a = 1 / (pow(da.x, 2) + pow(da.y, 2) + pow(da.z, 2));
            ra = sqrt(r2a);
            r6a = r2a * r2a * r2a;

            Vela = topo.coulomb_constant * qi * qj * ra;

            ai_aii = ai_type.aii_normal;
            aj_aii = aj_type.aii_normal;
            ai_bii = ai_type.bii_normal;
            aj_bii = aj_type.bii_normal;

            V_a = r6a * r6a * ai_aii * aj_aii;
            V_b = r6a * ai_bii * aj_bii;
            dva = r2a * (-Vela - 12 * V_a + 6 * V_b);

            dvelocities[i].x -= dva * da.x;
            dvelocities[i].y -= dva * da.y;
            dvelocities[i].z -= dva * da.z;

            dvelocities[j].x += dva * da.x;
            dvelocities[j].y += dva * da.y;
            dvelocities[j].z += dva * da.z;

            E_nonbond_pw.Ucoul += Vela;
            E_nonbond_pw.Uvdw += (V_a - V_b);

            // if (j - n_atoms_solute == 1754) {
            //     printf("MAT[%d][%d].W = %f %f %f\n", i, j - n_atoms_solute, dva *
            //     da.x, dva * da.y, dva * da.z);
            // }
        }
    }
}

void calc_nonbonded_pw_forces_host() {
    int mem_size_W = 3 * n_waters * sizeof(coord_t);
    int mem_size_DV_W = 3 * n_waters * sizeof(dvel_t);
    int mem_size_DV_X = n_atoms_solute * sizeof(dvel_t);
    int mem_size_PW_MAT = 3 * n_waters * n_patoms * sizeof(calc_pw_t);

    int n_blocks_p = (n_patoms + BLOCK_SIZE - 1) / BLOCK_SIZE;
    int n_blocks_w = (3 * n_waters + BLOCK_SIZE - 1) / BLOCK_SIZE;

    int mem_size_PW_Evdw = n_blocks_p * n_blocks_w * sizeof(double);
    int mem_size_PW_Ecoul = n_blocks_p * n_blocks_w * sizeof(double);

    if (!pw_gpu_set) {
        pw_gpu_set = true;

#ifdef DEBUG
    printf("Allocating PW_MAT\n");
#endif
        check_cudaMalloc((void **) &PW_MAT, mem_size_PW_MAT);

#ifdef DEBUG
    printf("Allocating PW_MAT\n");
#endif
        check_cudaMalloc((void **) &PW_MAT, mem_size_PW_MAT);
#ifdef DEBUG
    printf("Allocating PW_MAT\n");
#endif
        check_cudaMalloc((void **) &PW_MAT, mem_size_PW_MAT);

#ifdef DEBUG
    printf("Allocating D_PW_Evdw\n");
#endif
        check_cudaMalloc((void **) &D_PW_Evdw, mem_size_PW_Evdw);
#ifdef DEBUG
    printf("Allocating D_PW_Ecoul\n");
#endif
        check_cudaMalloc((void **) &D_PW_Ecoul, mem_size_PW_Ecoul);

        check_cudaMalloc((void **) &D_PW_evdw_TOT, sizeof(double));
        check_cudaMalloc((void **) &D_PW_ecoul_TOT, sizeof(double));

#ifdef DEBUG
    printf("All GPU solvent memory allocated\n");
#endif

        h_PW_MAT = (calc_pw_t *) malloc(mem_size_PW_MAT);
    }

    cudaMemcpy(W, &coords[n_atoms_solute], mem_size_W, cudaMemcpyHostToDevice);
    cudaMemcpy(DV_X, dvelocities, mem_size_DV_X, cudaMemcpyHostToDevice);
    cudaMemcpy(DV_W, &dvelocities[n_atoms_solute], mem_size_DV_W,
               cudaMemcpyHostToDevice);

    dim3 threads, grid;

    threads = dim3(BLOCK_SIZE, BLOCK_SIZE);
    grid = dim3((3 * n_waters + BLOCK_SIZE - 1) / threads.x,
                (n_patoms + BLOCK_SIZE - 1) / threads.y);

    calc_pw_dvel_matrix<<<grid, threads>>>(
        n_patoms, n_atoms_solute, n_waters, X, W, D_PW_Evdw, D_PW_Ecoul, PW_MAT,
        D_ccharges, D_charges, D_catypes, D_atypes, D_patoms, D_excluded, topo);
    calc_pw_dvel_vector_column<<<((3 * n_waters + BLOCK_SIZE - 1) / BLOCK_SIZE),
            BLOCK_SIZE>>>(n_patoms, n_waters, DV_X, DV_W,
                          PW_MAT);
    calc_pw_dvel_vector_row<<<((n_patoms + BLOCK_SIZE - 1) / BLOCK_SIZE),
            BLOCK_SIZE>>>(n_patoms, n_waters, DV_X, DV_W,
                          PW_MAT, D_patoms);

#ifdef DEBUG
  cudaMemcpy(h_PW_MAT, PW_MAT, mem_size_PW_MAT, cudaMemcpyDeviceToHost);
#endif

    // for (int i = 0; i < n_waters; i++) {
    //     printf("X[%d] = %f %f %f\n", i, coords[i].x, coords[i].y, coords[i].z);
    // }

    // printf("n_patoms = %d n_watoms = %d\n", n_patoms, 3 * n_waters);
#ifdef DEBUG
  for (int i = 0; i < n_patoms; i++) {
    for (int j = 0; j < 3 * n_waters; j++) {
      if (h_PW_MAT[3 * i * n_waters + j].P.x > 100)
        printf("PW_MAT[%d][%d].P = %f %f %f\n", i, j,
               h_PW_MAT[3 * i * n_waters + j].P.x,
               h_PW_MAT[3 * i * n_waters + j].P.y,
               h_PW_MAT[3 * i * n_waters + j].P.z);
    }
  }
#endif

    cudaMemcpy(dvelocities, DV_X, mem_size_DV_X, cudaMemcpyDeviceToHost);
    cudaMemcpy(&dvelocities[n_atoms_solute], DV_W, mem_size_DV_W,
               cudaMemcpyDeviceToHost);

    calc_energy_sum<<<1, threads>>>(n_blocks_p, n_blocks_w, D_PW_evdw_TOT,
                                    D_PW_ecoul_TOT, D_PW_Evdw, D_PW_Ecoul, false);

    cudaMemcpy(&PW_evdw_TOT, D_PW_evdw_TOT, sizeof(double),
               cudaMemcpyDeviceToHost);
    cudaMemcpy(&PW_ecoul_TOT, D_PW_ecoul_TOT, sizeof(double),
               cudaMemcpyDeviceToHost);

    E_nonbond_pw.Uvdw += PW_evdw_TOT;
    E_nonbond_pw.Ucoul += PW_ecoul_TOT;

    // for (int i = 0; i < n_atoms; i++) {
    //     printf("dvelocities[%d] = %f %f %f\n", i, dvelocities[i].x,
    //     dvelocities[i].y, dvelocities[i].z);
    // }
}

/* =============================================
 * == DEVICE
 * =============================================
 */

// W-W interactions

__device__ void set_water(int n_waters, int row, int column, dvel_t *val,
                          calc_ww_t *MAT) {
    MAT[column + n_waters * row].O = val[0];
    MAT[column + n_waters * row].H1 = val[1];
    MAT[column + n_waters * row].H2 = val[2];
}

__device__ void calc_ww_dvel_matrix_incr(int row, int column, double crg_ow,
                                         double crg_hw, double A_OO,
                                         double B_OO, coord_t *Xs, coord_t *Ys,
                                         double *Evdw, double *Ecoul,
                                         dvel_t *water_a, dvel_t *water_b,
                                         topo_t D_topo) {
    double rOX, rH1X, rH2X, r2;
    coord_t dOX, dH1X, dH2X;
    double Vel, V_a, V_b, dv;
    double tempX, tempY, tempZ;

    int i = 3 * row;
    int j = 3 * column;
    int wi = 0, wj = 0;

    // --- O - (O,H1,H2) ---
    dOX.x = Ys[j].x - Xs[i].x;
    dOX.y = Ys[j].y - Xs[i].y;
    dOX.z = Ys[j].z - Xs[i].z;
    rOX = pow(dOX.x, 2) + pow(dOX.y, 2) + pow(dOX.z, 2);

    // O - H1
    j += 1;
    wj += 1;

    dH1X.x = Ys[j].x - Xs[i].x;
    dH1X.y = Ys[j].y - Xs[i].y;
    dH1X.z = Ys[j].z - Xs[i].z;
    rH1X = pow(dH1X.x, 2) + pow(dH1X.y, 2) + pow(dH1X.z, 2);

    // O-H2 (X=H2)
    j += 1;
    wj += 1;

    dH2X.x = Ys[j].x - Xs[i].x;
    dH2X.y = Ys[j].y - Xs[i].y;
    dH2X.z = Ys[j].z - Xs[i].z;
    rH2X = pow(dH2X.x, 2) + pow(dH2X.y, 2) + pow(dH2X.z, 2);
    rOX = sqrt(1 / rOX);
    rH1X = sqrt(1 / rH1X);
    rH2X = sqrt(1 / rH2X);

    // O - O
    r2 = rOX * rOX;
    Vel = D_topo.coulomb_constant * pow(crg_ow, 2) * rOX;
    V_a = A_OO * (r2 * r2 * r2) * (r2 * r2 * r2);
    V_b = B_OO * (r2 * r2 * r2);
    *Evdw += (V_a - V_b);
    *Ecoul += Vel;
    dv = r2 * (-Vel - 12 * V_a + 6 * V_b);
    j -= 2; // move pointer back to O in interacting molecule
    wj -= 2;

    water_a[wi].x -= (dv * dOX.x);
    water_b[wj].x += (dv * dOX.x);
    water_a[wi].y -= (dv * dOX.y);
    water_b[wj].y += (dv * dOX.y);
    water_a[wi].z -= (dv * dOX.z);
    water_b[wj].z += (dv * dOX.z);

    // O - H1
    r2 = pow(rH1X, 2);
    Vel = D_topo.coulomb_constant * crg_ow * crg_hw * rH1X;
    *Ecoul += Vel;
    dv = r2 * (-Vel);
    j += 1; // point to H1 in j-molecule
    wj += 1;

    water_a[wi].x -= (dv * dH1X.x);
    water_b[wj].x += (dv * dH1X.x);
    water_a[wi].y -= (dv * dH1X.y);
    water_b[wj].y += (dv * dH1X.y);
    water_a[wi].z -= (dv * dH1X.z);
    water_b[wj].z += (dv * dH1X.z);

    // O - H2
    r2 = pow(rH2X, 2);
    Vel = D_topo.coulomb_constant * crg_ow * crg_hw * rH2X;
    *Ecoul += Vel;
    dv = r2 * (-Vel);
    j += 1; // point to H2 in j-molecule
    wj += 1;

    water_a[wi].x -= (dv * dH2X.x);
    water_b[wj].x += (dv * dH2X.x);
    water_a[wi].y -= (dv * dH2X.y);
    water_b[wj].y += (dv * dH2X.y);
    water_a[wi].z -= (dv * dH2X.z);
    water_b[wj].z += (dv * dH2X.z);

    // --- H1 - (O,H1,H2) ---
    i += 1; // Point to H1 in i-molecule
    wi += 1;
    j -= 2; // Point to O in j-molecule
    wj -= 2;

    // H1 - O (X=O)
    dOX.x = Ys[j].x - Xs[i].x;
    dOX.y = Ys[j].y - Xs[i].y;
    dOX.z = Ys[j].z - Xs[i].z;
    rOX = pow(dOX.x, 2) + pow(dOX.y, 2) + pow(dOX.z, 2);

    // H1 - H1 (X=H1)
    j += 1;
    wj += 1;

    dH1X.x = Ys[j].x - Xs[i].x;
    dH1X.y = Ys[j].y - Xs[i].y;
    dH1X.z = Ys[j].z - Xs[i].z;
    rH1X = pow(dH1X.x, 2) + pow(dH1X.y, 2) + pow(dH1X.z, 2);

    // H1 - H2 (X=H2)
    j += 1;
    wj += 1;

    dH2X.x = Ys[j].x - Xs[i].x;
    dH2X.y = Ys[j].y - Xs[i].y;
    dH2X.z = Ys[j].z - Xs[i].z;
    rH2X = pow(dH2X.x, 2) + pow(dH2X.y, 2) + pow(dH2X.z, 2);
    rOX = sqrt(1 / rOX);
    rH1X = sqrt(1 / rH1X);
    rH2X = sqrt(1 / rH2X);

    // H1 - O
    r2 = rOX * rOX;
    Vel = D_topo.coulomb_constant * crg_hw * crg_ow * rOX;
    *Ecoul += Vel;
    dv = r2 * (-Vel);
    j -= 2; // move pointer back to O in interacting molecule
    wj -= 2;
    water_a[wi].x -= (dv * dOX.x);
    water_b[wj].x += (dv * dOX.x);
    water_a[wi].y -= (dv * dOX.y);
    water_b[wj].y += (dv * dOX.y);
    water_a[wi].z -= (dv * dOX.z);
    water_b[wj].z += (dv * dOX.z);

    // H1 - H1
    r2 = pow(rH1X, 2);
    Vel = D_topo.coulomb_constant * crg_hw * crg_hw * rH1X;
    *Ecoul += Vel;
    dv = r2 * (-Vel);
    j += 1; // point to H1 in j-molecule
    wj += 1;
    water_a[wi].x -= (dv * dH1X.x);
    water_b[wj].x += (dv * dH1X.x);
    water_a[wi].y -= (dv * dH1X.y);
    water_b[wj].y += (dv * dH1X.y);
    water_a[wi].z -= (dv * dH1X.z);
    water_b[wj].z += (dv * dH1X.z);

    // H1 - H2
    r2 = pow(rH2X, 2);
    Vel = D_topo.coulomb_constant * crg_hw * crg_hw * rH2X;
    *Ecoul += Vel;
    dv = r2 * (-Vel);
    j += 1; // point to H2 in j-molecule
    wj += 1;
    water_a[wi].x -= (dv * dH2X.x);
    water_b[wj].x += (dv * dH2X.x);
    water_a[wi].y -= (dv * dH2X.y);
    water_b[wj].y += (dv * dH2X.y);
    water_a[wi].z -= (dv * dH2X.z);
    water_b[wj].z += (dv * dH2X.z);

    // --- H2 - (O,H1,H2) ---
    i += 1; // Point to H2 in i-molecule
    wi += 1;
    j -= 2; // Point to O in j-molecule
    wj -= 2;

    // H2 - O (X=O)
    dOX.x = Ys[j].x - Xs[i].x;
    dOX.y = Ys[j].y - Xs[i].y;
    dOX.z = Ys[j].z - Xs[i].z;
    rOX = pow(dOX.x, 2) + pow(dOX.y, 2) + pow(dOX.z, 2);

    // H2 - H1 (X=H1)
    j += 1;
    wj += 1;

    dH1X.x = Ys[j].x - Xs[i].x;
    dH1X.y = Ys[j].y - Xs[i].y;
    dH1X.z = Ys[j].z - Xs[i].z;
    rH1X = pow(dH1X.x, 2) + pow(dH1X.y, 2) + pow(dH1X.z, 2);

    // H2 - H2 (X=H2)
    j += 1;
    wj += 1;

    dH2X.x = Ys[j].x - Xs[i].x;
    dH2X.y = Ys[j].y - Xs[i].y;
    dH2X.z = Ys[j].z - Xs[i].z;
    rH2X = pow(dH2X.x, 2) + pow(dH2X.y, 2) + pow(dH2X.z, 2);
    rOX = sqrt(1 / rOX);
    rH1X = sqrt(1 / rH1X);
    rH2X = sqrt(1 / rH2X);

    // H2 - O
    r2 = rOX * rOX;
    Vel = D_topo.coulomb_constant * crg_hw * crg_ow * rOX;
    *Ecoul += Vel;
    dv = r2 * (-Vel);
    j -= 2; // move pointer back to O in interacting molecule
    wj -= 2;
    water_a[wi].x -= (dv * dOX.x);
    water_b[wj].x += (dv * dOX.x);
    water_a[wi].y -= (dv * dOX.y);
    water_b[wj].y += (dv * dOX.y);
    water_a[wi].z -= (dv * dOX.z);
    water_b[wj].z += (dv * dOX.z);

    // H2 - H1
    r2 = pow(rH1X, 2);
    Vel = D_topo.coulomb_constant * crg_hw * crg_hw * rH1X;
    *Ecoul += Vel;
    dv = r2 * (-Vel);
    j += 1; // point to H1 in j-molecule
    wj += 1;
    water_a[wi].x -= (dv * dH1X.x);
    water_b[wj].x += (dv * dH1X.x);
    water_a[wi].y -= (dv * dH1X.y);
    water_b[wj].y += (dv * dH1X.y);
    water_a[wi].z -= (dv * dH1X.z);
    water_b[wj].z += (dv * dH1X.z);

    // H1 - H2
    r2 = pow(rH2X, 2);
    Vel = D_topo.coulomb_constant * crg_hw * crg_hw * rH2X;
    *Ecoul += Vel;
    dv = r2 * (-Vel);
    j += 1; // point to H2 in j-molecule
    wj += 1;
    water_a[wi].x -= (dv * dH2X.x);
    water_b[wj].x += (dv * dH2X.x);
    water_a[wi].y -= (dv * dH2X.y);
    water_b[wj].y += (dv * dH2X.y);
    water_a[wi].z -= (dv * dH2X.z);
    water_b[wj].z += (dv * dH2X.z);
}

__global__ void calc_ww_dvel_matrix(int n_waters, double crg_ow, double crg_hw,
                                    double A_OO, double B_OO, coord_t *W,
                                    double *Evdw, double *Ecoul, calc_ww_t *MAT,
                                    topo_t D_topo) {
    // Block index
    int bx = blockIdx.x;
    int by = blockIdx.y;

    // Thread index
    int tx = threadIdx.x;
    int ty = threadIdx.y;

    __shared__ double Ecoul_S[BLOCK_SIZE][BLOCK_SIZE];
    __shared__ double Evdw_S[BLOCK_SIZE][BLOCK_SIZE];

    Ecoul_S[ty][tx] = 0;
    Evdw_S[ty][tx] = 0;

    // if (bx == 9 && by == 0) printf("bx = %d by = %d tx = %d ty = %d\n", bx, by,
    // tx, ty);

    int aStart = 3 * BLOCK_SIZE * by;
    int bStart = 3 * BLOCK_SIZE * bx;

    if (aStart + 3 * ty >= 3 * n_waters)
        return;
    if (bStart + 3 * tx >= 3 * n_waters)
        return;

    __shared__ coord_t Xs[3 * BLOCK_SIZE];
    __shared__ coord_t Ys[3 * BLOCK_SIZE];

    if (tx == 0) {
        Xs[3 * ty] = W[aStart + 3 * ty];
        Xs[3 * ty + 1] = W[aStart + 3 * ty + 1];
        Xs[3 * ty + 2] = W[aStart + 3 * ty + 2];
    }

    if (ty == 0) {
        Ys[3 * tx] = W[bStart + 3 * tx];
        Ys[3 * tx + 1] = W[bStart + 3 * tx + 1];
        Ys[3 * tx + 2] = W[bStart + 3 * tx + 2];
    }

    __syncthreads();

    if (bx < by || (bx == by && tx < ty))
        return;

    dvel_t water_a[3], water_b[3];
    memset(&water_a, 0, 3 * sizeof(dvel_t));
    memset(&water_b, 0, 3 * sizeof(dvel_t));

    int row = by * BLOCK_SIZE + ty;
    int column = bx * BLOCK_SIZE + tx;

    // if (row == 0 && column == 1) {
    //     printf("Xs[0] = %f\n", Xs[0]);
    //     printf("Ys[0] = %f\n", Ys[0]);
    //     printf("Xs[1] = %f\n", Xs[1]);
    //     printf("Ys[1] = %f\n", Ys[1]);
    //     printf("Xs[2] = %f\n", Xs[2]);
    //     printf("Ys[2] = %f\n", Ys[2]);
    //     printf("Xs[3] = %f\n", Xs[3]);
    //     printf("Ys[3] = %f\n", Ys[3]);
    //     printf("Xs[4] = %f\n", Xs[4]);
    //     printf("Ys[4] = %f\n", Ys[4]);

    //     printf("Ys[%d] = %f Xs[%d] = %f\n", 3 * ty, Ys[3 * ty], 3 * tx, Xs[3 *
    //     tx]);
    // }

    if (bx != by || tx != ty) {
        double evdw = 0, ecoul = 0;
        calc_ww_dvel_matrix_incr(ty, tx, crg_ow, crg_hw, A_OO, B_OO, Xs, Ys, &evdw,
                                 &ecoul, water_a, water_b, D_topo);
        Evdw_S[ty][tx] = evdw;
        Ecoul_S[ty][tx] = ecoul;
        // if (bx == 2 && tx == 0 && by == 1 && ty == 0) {
        //     printf("evdw = %f, ecoul = %f, \n", evdw, ecoul);
        // }
    }

    // if (row == 0 && column == 1) {
    //     printf("water_a = %f %f %f water_b = %f %f %f\n", water_a[0].x,
    //     water_a[0].y, water_a[0].z, water_b[0].x, water_b[0].y, water_b[0].z);
    // }

    set_water(n_waters, row, column, water_a, MAT);
    set_water(n_waters, column, row, water_b, MAT);

    __syncthreads();

    int rowlen = (n_waters + BLOCK_SIZE - 1) / BLOCK_SIZE;

    if (tx == 0 && ty == 0) {
        double tot_Evdw = 0;
        double tot_Ecoul = 0;
        for (int i = 0; i < BLOCK_SIZE; i++) {
            for (int j = 0; j < BLOCK_SIZE; j++) {
                tot_Evdw += Evdw_S[i][j];
                tot_Ecoul += Ecoul_S[i][j];
            }
        }
        Evdw[rowlen * by + bx] = tot_Evdw;
        Ecoul[rowlen * by + bx] = tot_Ecoul;

        // if (bx == 70 && by == 70) {
        //     printf("bx = %d, by = %d, tot_Evdw = %f, tot_Ecoul = %f\n", bx, by,
        //     tot_Evdw, tot_Ecoul);
        // }
        // printf("bx = %d by = %d tot_Evdw = %f tot_Ecoul = %f\n", bx, by,
        // Evdw[rowlen * by + bx], Ecoul[rowlen * by + bx]);
    }

    // if (bx == 8 && by == 2 && tx == 0 && ty == 0) {
    //     for (int i = 0; i < BLOCK_SIZE; i++) {
    //         for (int j = 0; j < BLOCK_SIZE; j++) {
    //             printf("Evdw_S[%d][%d] = %f\n", i, j, Evdw_S[i][j]);
    //             printf("Ecoul_S[%d][%d] = %f\n", i, j, Ecoul_S[i][j]);
    //         }
    //     }
    //     printf("tot_Evdw = %f\n", Evdw[rowlen * by + bx]);
    //     printf("tot_Ecoul = %f\n", Ecoul[rowlen * by + bx]);
    // }
    __syncthreads();
}

__global__ void calc_ww_dvel_vector(int n_waters, dvel_t *DV_W,
                                    calc_ww_t *MAT) {
    int row = blockIdx.x * blockDim.x + threadIdx.x;
    if (row >= n_waters)
        return;

    dvel_t dO, dH1, dH2;
    dO.x = 0;
    dO.y = 0;
    dO.z = 0;
    dH1.x = 0;
    dH1.y = 0;
    dH1.z = 0;
    dH2.x = 0;
    dH2.y = 0;
    dH2.z = 0;

    for (int i = 0; i < n_waters; i++) {
        if (i != row) {
            dO.x += MAT[i + n_waters * row].O.x;
            dO.y += MAT[i + n_waters * row].O.y;
            dO.z += MAT[i + n_waters * row].O.z;
            dH1.x += MAT[i + n_waters * row].H1.x;
            dH1.y += MAT[i + n_waters * row].H1.y;
            dH1.z += MAT[i + n_waters * row].H1.z;
            dH2.x += MAT[i + n_waters * row].H2.x;
            dH2.y += MAT[i + n_waters * row].H2.y;
            dH2.z += MAT[i + n_waters * row].H2.z;
        }
    }

    DV_W[3 * row].x += dO.x;
    DV_W[3 * row].y += dO.y;
    DV_W[3 * row].z += dO.z;
    DV_W[3 * row + 1].x += dH1.x;
    DV_W[3 * row + 1].y += dH1.y;
    DV_W[3 * row + 1].z += dH1.z;
    DV_W[3 * row + 2].x += dH2.x;
    DV_W[3 * row + 2].y += dH2.y;
    DV_W[3 * row + 2].z += dH2.z;

    __syncthreads();
}

// P-W interactions

__device__ void calc_pw_dvel_matrix_incr(
    int row, int pi, int column, int j, int n_atoms_solute, coord_t *Xs,
    coord_t *Ws, bool *excluded_s, double *Evdw, double *Ecoul, calc_pw_t *pw,
    ccharge_t *D_ccharges, charge_t *D_charges, catype_t *D_catypes,
    atype_t *D_atypes, p_atom_t *D_patoms, topo_t D_topo) {
    coord_t da;
    double r2a, ra, r6a;
    double Vela, V_a, V_b;
    double dva;
    double qi, qj;
    double ai_aii, aj_aii, ai_bii, aj_bii;
    catype_t ai_type, aj_type;
    atype_t i_type, j_type;

    if (excluded_s[row])
        return;
    qi = D_ccharges[D_charges[pi].code - 1].charge;
    qj = D_ccharges[D_charges[n_atoms_solute + j].code - 1]
            .charge; // TODO: FIX THIS!!! WILL NOT WORK WITH QATOMS!!!!!

    // if (pi < 100 && j < 100){
    //     printf("qi = %f qj = %f\n", qi, qj);
    // }

    ai_type = D_catypes[D_atypes[pi].code - 1];
    aj_type = D_catypes[D_atypes[n_atoms_solute + j].code - 1];

    da.x = Ws[column].x - Xs[row].x;
    da.y = Ws[column].y - Xs[row].y;
    da.z = Ws[column].z - Xs[row].z;
    r2a = 1 / (pow(da.x, 2) + pow(da.y, 2) + pow(da.z, 2));
    ra = sqrt(r2a);
    r6a = r2a * r2a * r2a;

    Vela = D_topo.coulomb_constant * qi * qj * ra;

    ai_aii = ai_type.aii_normal;
    aj_aii = aj_type.aii_normal;
    ai_bii = ai_type.bii_normal;
    aj_bii = aj_type.bii_normal;

    V_a = r6a * r6a * ai_aii * aj_aii;
    V_b = r6a * ai_bii * aj_bii;
    dva = r2a * (-Vela - 12 * V_a + 6 * V_b);

    pw->P.x -= dva * da.x;
    pw->P.y -= dva * da.y;
    pw->P.z -= dva * da.z;

    pw->W.x += dva * da.x;
    pw->W.y += dva * da.y;
    pw->W.z += dva * da.z;

    *Ecoul += Vela;
    *Evdw += (V_a - V_b);

    // if (pi == 522 && j == 175) {
    //     printf("Vela = %f V_a = %f V_b = %f P = %f %f %f ai_aii = %f aj_aii =
    //     %f\n", Vela, V_a, V_b, pw->P.x, pw->P.y, pw->P.z, ai_aii, aj_aii);
    // }

    // if (pi < 100 && j < 100) printf("Evdw = %f Ecoul = %f\n", *Evdw, *Ecoul);
}

__global__ void calc_pw_dvel_matrix(int n_patoms, int n_atoms_solute,
                                    int n_waters, coord_t *X, coord_t *W,
                                    double *Evdw, double *Ecoul,
                                    calc_pw_t *PW_MAT, ccharge_t *D_ccharges,
                                    charge_t *D_charges, catype_t *D_catypes,
                                    atype_t *D_atypes, p_atom_t *D_patoms,
                                    bool *D_excluded, topo_t D_topo) {
    // Block index
    int bx = blockIdx.x;
    int by = blockIdx.y;

    // Thread index
    int tx = threadIdx.x;
    int ty = threadIdx.y;

    __shared__ double Ecoul_S[BLOCK_SIZE][BLOCK_SIZE];
    __shared__ double Evdw_S[BLOCK_SIZE][BLOCK_SIZE];

    Ecoul_S[ty][tx] = 0;
    Evdw_S[ty][tx] = 0;

    // if (tx == 0 && ty == 0) printf("bx = %d by = %d\n", bx, by);

    // if (bx == 0 && by == 0) printf("bx = %d by = %d tx = %d ty = %d\n", bx, by,
    // tx, ty);

    int aStart = BLOCK_SIZE * by;
    int bStart = BLOCK_SIZE * bx;

    if (aStart + ty >= n_patoms)
        return;
    if (bStart + tx >= 3 * n_waters)
        return;

    // if (bx == 8 && by == 1) printf("bx = %d by = %d tx = %d ty = %d\n", bx, by,
    // tx, ty);

    __shared__ coord_t Xs[BLOCK_SIZE];
    __shared__ coord_t Ws[BLOCK_SIZE];
    __shared__ bool excluded_s[BLOCK_SIZE];

    int pi = D_patoms[aStart + ty].a - 1;

    Xs[ty] = X[pi];
    Ws[tx] = W[bStart + tx];

    if (tx == 0) {
        excluded_s[ty] = D_excluded[pi];
    }

    __syncthreads();

    calc_pw_t pw;
    memset(&pw, 0, sizeof(calc_pw_t));

    int row = by * BLOCK_SIZE + ty;
    int column = bx * BLOCK_SIZE + tx;

    // if (row == 0 && column == 1) {
    //     printf("Xs[0] = %f\n", Xs[0]);
    //     printf("Ys[0] = %f\n", Ys[0]);
    //     printf("Xs[1] = %f\n", Xs[1]);
    //     printf("Ys[1] = %f\n", Ys[1]);
    //     printf("Xs[2] = %f\n", Xs[2]);
    //     printf("Ys[2] = %f\n", Ys[2]);
    //     printf("Xs[3] = %f\n", Xs[3]);
    //     printf("Ys[3] = %f\n", Ys[3]);
    //     printf("Xs[4] = %f\n", Xs[4]);
    //     printf("Ys[4] = %f\n", Ys[4]);

    //     printf("Ys[%d] = %f Xs[%d] = %f\n", 3 * ty, Ys[3 * ty], 3 * tx, Xs[3 *
    //     tx]);
    // }

    // if (bx == 8 && by == 1) printf("bx = %d by = %d tx = %d ty = %d\n", bx, by,
    // tx, ty);
    // __device__ void calc_pw_dvel_matrix_incr(int row, int pi, int column, int
    // j, int n_patoms, coord_t *Ps, coord_t *Xs, double *Evdw, double *Ecoul,
    // calc_pw_t *pw, ccharge_t *D_ccharges, charge_t *D_charges, catype_t
    // *D_catypes, atype_t *D_atypes, p_atom_t *D_patoms)
    double evdw = 0, ecoul = 0;
    calc_pw_dvel_matrix_incr(ty, pi, tx, bStart + tx, n_atoms_solute, Xs, Ws,
                             excluded_s, &evdw, &ecoul, &pw, D_ccharges,
                             D_charges, D_catypes, D_atypes, D_patoms, D_topo);
    Evdw_S[ty][tx] = evdw;
    Ecoul_S[ty][tx] = ecoul;

    // if (row == 0 && column == 1) {
    //     printf("water_a = %f %f %f water_b = %f %f %f\n", water_a[0].x,
    //     water_a[0].y, water_a[0].z, water_b[0].x, water_b[0].y, water_b[0].z);
    // }

    // if (bx == 8 && by == 1) printf("n_qatoms = %d\n", n_qatoms);
    // if (bx == 8 && by == 1) printf("qi = %d j = %d charge[%d] = %f\n", row,
    // column, row + n_qatoms, D_qcharges[row + n_qatoms * 1].q);

    PW_MAT[column + 3 * n_waters * row] = pw;

    __syncthreads();

    int rowlen = (3 * n_waters + BLOCK_SIZE - 1) / BLOCK_SIZE;

    if (tx == 0 && ty == 0) {
        double tot_Evdw = 0;
        double tot_Ecoul = 0;
        for (int i = 0; i < BLOCK_SIZE; i++) {
            for (int j = 0; j < BLOCK_SIZE; j++) {
                tot_Evdw += Evdw_S[i][j];
                tot_Ecoul += Ecoul_S[i][j];
            }
        }
        Evdw[rowlen * by + bx] = tot_Evdw;
        Ecoul[rowlen * by + bx] = tot_Ecoul;
    }

    __syncthreads();
}

__global__ void calc_pw_dvel_vector_row(int n_patoms, int n_waters,
                                        dvel_t *DV_X, dvel_t *DV_W,
                                        calc_pw_t *PW_MAT, p_atom_t *D_patoms) {
    int row = blockIdx.x * blockDim.x + threadIdx.x;
    if (row >= n_patoms)
        return;

    dvel_t dP;

    dP.x = 0;
    dP.y = 0;
    dP.z = 0;

    for (int i = 0; i < 3 * n_waters; i++) {
        dP.x += PW_MAT[i + 3 * n_waters * row].P.x;
        dP.y += PW_MAT[i + 3 * n_waters * row].P.y;
        dP.z += PW_MAT[i + 3 * n_waters * row].P.z;
    }

    int p = D_patoms[row].a - 1;

    DV_X[p].x += dP.x;
    DV_X[p].y += dP.y;
    DV_X[p].z += dP.z;

    __syncthreads();
}

__global__ void calc_pw_dvel_vector_column(int n_patoms, int n_waters,
                                           dvel_t *DV_X, dvel_t *DV_W,
                                           calc_pw_t *PW_MAT) {
    int column = blockIdx.x * blockDim.x + threadIdx.x;
    if (column >= 3 * n_waters)
        return;

    dvel_t dW;

    dW.x = 0;
    dW.y = 0;
    dW.z = 0;

    for (int i = 0; i < n_patoms; i++) {
        dW.x += PW_MAT[column + 3 * n_waters * i].W.x;
        dW.y += PW_MAT[column + 3 * n_waters * i].W.y;
        dW.z += PW_MAT[column + 3 * n_waters * i].W.z;
    }

    DV_W[column].x += dW.x;
    DV_W[column].y += dW.y;
    DV_W[column].z += dW.z;

    __syncthreads();
}

// General
__global__ void calc_energy_sum(int rows, int columns, double *Evdw_TOT,
                                double *Ecoul_TOT, double *Evdw, double *Ecoul,
                                bool upper_diagonal) {
    int tx = threadIdx.x;
    int ty = threadIdx.y;

    __shared__ double Ecoul_S[BLOCK_SIZE][BLOCK_SIZE];
    __shared__ double Evdw_S[BLOCK_SIZE][BLOCK_SIZE];

    // TODO: better way to distribute upper diagonal over threads? Seems like
    // threads in left bottom corner have less work
    double coul_TOT = 0;
    double vdw_TOT = 0;
    for (int i = ty; i < rows; i += BLOCK_SIZE) {
        for (int j = tx; j < columns; j += BLOCK_SIZE) {
            if (i <= j || !upper_diagonal) {
                coul_TOT += Ecoul[i * columns + j];
                vdw_TOT += Evdw[i * columns + j];
            }
        }
    }
    Ecoul_S[ty][tx] = coul_TOT;
    Evdw_S[ty][tx] = vdw_TOT;

    __syncthreads();

    if (tx == 0 && ty == 0) {
        double Evdw_temp = 0;
        double Ecoul_temp = 0;

        for (int i = 0; i < BLOCK_SIZE; i++) {
            for (int j = 0; j < BLOCK_SIZE; j++) {
                Evdw_temp += Evdw_S[i][j];
                Ecoul_temp += Ecoul_S[i][j];
            }
        }

        *Evdw_TOT = Evdw_temp;
        *Ecoul_TOT = Ecoul_temp;
    }
}

void clean_d_solvent() {
    cudaFree(W);
    cudaFree(WW_MAT);
    cudaFree(PW_MAT);
    cudaFree(DV_W);

    free(h_WW_MAT);
    free(h_PW_MAT);

#ifdef DEBUG
  printf("All memory freed\n");
#endif
}