#include "cuda/include/CudaContext.cuh"
#include "cuda/include/CudaNonbondedWWForce.cuh"
#include <iostream>
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

__global__ void calc_ww(const int N, const double crg_ow, const double crg_hw,
                        const double A_OO, const double B_OO, const topo_t D_topo,
                        coord_t* __restrict__ W, dvel_t* __restrict__ DV_W,
                        double* __restrict__ Evdw_TOT, double* __restrict__ ecoul_TOT) {
    const int warpsPerBlock = blockDim.x >> 5;  // THREAD_BLOCK_SIZE / 32
    const int tid = threadIdx.x;
    const int lane = tid & 31;
    const int warpInBlock = tid >> 5;

    const int tile = blockIdx.x * warpsPerBlock + warpInBlock;

    // Map tile -> (x,y) with y>=x (upper triangle)
    int x, y;
    idx2xy((N + 31) >> 5, tile, x, y);  // nBlocks = ceil(N/32)

    const int baseX = x << 5;  // x*32
    const int baseY = y << 5;  // y*32

    const int rowIdx = baseY + lane;
    const int colIdx0 = baseX + lane;  // lane's "original atom2"

    coord_t invalid = {-1e9, -1e9, -1e9};
    coord_t row = (rowIdx < N ? W[rowIdx] : invalid);

    // "col" state that will be rotated around the warp
    int colIdx = colIdx0;
    coord_t col = (colIdx < N ? W[colIdx] : invalid);

    // Accumulators
    double3 rowForce = {0.0, 0.0, 0.0};
    double3 colForce = {0.0, 0.0, 0.0};

    double evdw_sum = 0.0;
    double ecoul_sum = 0.0;

    const unsigned mask = 0xffffffffu;

// 32-step warp rotation = warp-synchronous "reduction" for atom2
#pragma unroll
    for (int i = 0; i < 32; i++) {
        // Skip invalid lanes / out-of-range atoms
        if (rowIdx < N && colIdx < N) {
            // For diagonal tile, only compute each pair once
            if (x != y || colIdx > rowIdx) {
                double evdw = 0.0, ecoul = 0.0, dv = 0.0, tx, ty, tz;
                calculate_unforce_bound(rowIdx, colIdx, row, col, D_topo,
                                        crg_ow, crg_hw, A_OO, B_OO,
                                        evdw, ecoul, dv, tx, ty, tz);

                evdw_sum += evdw;
                ecoul_sum += ecoul;
                double v_x = dv * tx;
                double v_y = dv * ty;
                double v_z = dv * tz;

                colForce.x += v_x;
                colForce.y += v_y;
                colForce.z += v_z;

                rowForce.x -= v_x;
                rowForce.y -= v_y;
                rowForce.z -= v_z;
            }
        }

        // Rotate (colIdx, col, colForce) so the "atom2 packet" visits every lane.
        const int src = (lane + 1) & 31;  // cyclic shift
        colIdx = __shfl_sync(mask, colIdx, src);
        col = shfl_coord(col, src, mask);
        colForce.x = shfl(colForce.x, src, mask);
        colForce.y = shfl(colForce.y, src, mask);
        colForce.z = shfl(colForce.z, src, mask);
    }

    // Write forces: one atomic per atom (row + "original atom2")
    if (rowIdx < N) {
        atomicAdd(&DV_W[rowIdx].x, rowForce.x);
        atomicAdd(&DV_W[rowIdx].y, rowForce.y);
        atomicAdd(&DV_W[rowIdx].z, rowForce.z);
    }
    if (colIdx0 < N) {
        atomicAdd(&DV_W[colIdx0].x, colForce.x);
        atomicAdd(&DV_W[colIdx0].y, colForce.y);
        atomicAdd(&DV_W[colIdx0].z, colForce.z);
    }

    for (int offset = 16; offset > 0; offset >>= 1) {
        evdw_sum += __shfl_down_sync(0xffffffffu, evdw_sum, offset);
        ecoul_sum += __shfl_down_sync(0xffffffffu, ecoul_sum, offset);
    }
    // printf("evdw_sum: %f, ecoul_sum: %f\n", evdw_sum, ecoul_sum);
    if (lane == 0) {
        atomicAdd(Evdw_TOT, evdw_sum);
        atomicAdd(ecoul_TOT, ecoul_sum);
    }
}
}  // namespace CudaNonbondedWWForce
void calc_nonbonded_ww_forces_host_v2() {
    using namespace CudaNonbondedWWForce;
    int N = 3 * n_waters;
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

    const int thread_num = 128;
    dim3 block_sz = dim3(thread_num);

    int tile_num_per_block = thread_num / 32;
    int nBlocks = (N + 31) / 32;
    int total_tiles = nBlocks * (nBlocks + 1) / 2;

    int grid_sz = (total_tiles + tile_num_per_block - 1) / tile_num_per_block;
    dim3 grid = dim3(grid_sz);

    calc_ww<<<grid, block_sz>>>(N, crg_ow, crg_hw, A_OO, B_OO, topo, W, DV_W,
                                D_WW_evdw_TOT, D_WW_ecoul_TOT);

    cudaDeviceSynchronize();
    cudaMemcpy(&dvelocities[n_atoms_solute], DV_W, mem_size_DV_W,
               cudaMemcpyDeviceToHost);
    cudaMemcpy(&WW_evdw_TOT, D_WW_evdw_TOT, sizeof(double),
               cudaMemcpyDeviceToHost);
    cudaMemcpy(&WW_ecoul_TOT, D_WW_ecoul_TOT, sizeof(double),
               cudaMemcpyDeviceToHost);
    
    


    // printf("WW E_vdw: %f, E_coul: %f\n", WW_evdw_TOT, WW_ecoul_TOT);
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