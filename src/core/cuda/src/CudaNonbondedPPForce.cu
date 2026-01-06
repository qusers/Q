#include <iostream>

#include "cuda/include/CudaContext.cuh"
#include "cuda/include/CudaNonbondedPPForce.cuh"
#include "utils.h"

namespace CudaNonbondedPPForce {
bool is_initialized = false;
double *D_PP_evdw_TOT, *D_PP_ecoul_TOT;

__device__ void calc_pp_dvel_matrix_incr(int row, int pi, int column, int pj,
                                         coord_t* Xs, coord_t* Ys, int* LJs, bool* excluded_s, double* Evdw, double* Ecoul, dvel_t* patom_a, dvel_t* patom_b,
                                         ccharge_t* D_ccharges, charge_t* D_charges, catype_t* D_catypes, atype_t* D_atypes, p_atom_t* D_patoms, topo_t D_topo) {
    bool bond14, bond23;
    double scaling;
    coord_t da;
    double r2a, ra, r6a;
    double Vela, V_a, V_b;
    double dva;
    double crg_i, crg_j;
    double ai_aii, aj_aii, ai_bii, aj_bii;
    catype_t ai_type, aj_type;

    bond23 = LJs[row * BLOCK_SIZE + column] == 3;
    bond14 = LJs[row * BLOCK_SIZE + column] == 1;

    if (bond23) return;
    if (excluded_s[row] || excluded_s[BLOCK_SIZE + column]) return;

    scaling = bond14 ? D_topo.el14_scale : 1;

    crg_i = D_ccharges[D_charges[pi].code - 1].charge;
    crg_j = D_ccharges[D_charges[pj].code - 1].charge;

    ai_type = D_catypes[D_atypes[pi].code - 1];
    aj_type = D_catypes[D_atypes[pj].code - 1];

    da.x = Ys[column].x - Xs[row].x;
    da.y = Ys[column].y - Xs[row].y;
    da.z = Ys[column].z - Xs[row].z;
    r2a = 1 / (pow(da.x, 2) + pow(da.y, 2) + pow(da.z, 2));
    ra = sqrt(r2a);
    r6a = r2a * r2a * r2a;

    Vela = scaling * D_topo.coulomb_constant * crg_i * crg_j * ra;

    ai_aii = bond14 ? ai_type.aii_1_4 : ai_type.aii_normal;
    aj_aii = bond14 ? aj_type.aii_1_4 : aj_type.aii_normal;
    ai_bii = bond14 ? ai_type.bii_1_4 : ai_type.bii_normal;
    aj_bii = bond14 ? aj_type.bii_1_4 : aj_type.bii_normal;

    V_a = r6a * r6a * ai_aii * aj_aii;
    V_b = r6a * ai_bii * aj_bii;
    dva = r2a * (-Vela - 12 * V_a + 6 * V_b);

    patom_a->x = -dva * da.x;
    patom_a->y = -dva * da.y;
    patom_a->z = -dva * da.z;

    patom_b->x = dva * da.x;
    patom_b->y = dva * da.y;
    patom_b->z = dva * da.z;

    *Ecoul += Vela;
    *Evdw += (V_a - V_b);
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

__device__ void calculate_unforce_bound(
    // const int x_idx,
    // const int y_idx,

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

    double& evdw,
    double& ecoul,
    double& dv) {
    double3 d = {y.x - x.x, y.y - x.y, y.z - x.z};

    double r2 = 1.0 / (d.x * d.x + d.y * d.y + d.z * d.z);
    double r = sqrt(r2);
    double r6 = r2 * r2 * r2;

    ecoul = scaling * coulomb_constant * x_charge * y_charge * r;

    double v_a = r6 * r6 * x_aii * y_aii;
    double v_b = r6 * x_bii * y_bii;
    evdw = v_a - v_b;
    dv = r2 * (-ecoul - 12.0 * v_a + 6.0 * v_b);
}

__global__ void calc_pp(
    const int N,
    const int n_atoms_solute,

    const charge_t* D_charges,    // charge type of each atom
    const ccharge_t* D_ccharges,  // charge value of each charge type

    const atype_t* D_atypes,
    const catype_t* D_catypes,

    const topo_t D_topo,

    const bool* D_excluded,

    const int* D_LJ_matrix,

    const p_atom_t* D_patoms,

    coord_t* X,
    dvel_t* DV_X,

    double* Evdw_TOT,
    double* Ecoul_TOT) {
    const int nBlocks = (N + 31) >> 5;  // ceil(N/32)
    const int totalTiles = nBlocks * (nBlocks + 1) / 2;

    const int warpsPerBlock = blockDim.x >> 5;  // THREAD_BLOCK_SIZE / 32
    const int tid = threadIdx.x;
    const int lane = tid & 31;
    const int warpInBlock = tid >> 5;

    const int tile = blockIdx.x * warpsPerBlock + warpInBlock;
    if (tile >= totalTiles) {
        return;  // extra warps in last block
    }

    // Map tile -> (x,y) with y>=x (upper triangle)
    int tile_x, tile_y;
    idx2xy(nBlocks, tile, tile_x, tile_y);

    const int baseX = tile_x << 5;  // x*32
    const int baseY = tile_y << 5;  // y*32

    int x_idx = baseX + lane;
    int y_idx = baseY + lane;

    int x_atom = (x_idx < N) ? (D_patoms[x_idx].a - 1) : -1;
    int y_atom = (y_idx < N) ? (D_patoms[y_idx].a - 1) : -1;

    coord_t invalid = {-1e9, -1e9, -1e9};
    coord_t x_coord = (x_atom >= 0) ? X[x_atom] : invalid;
    coord_t y_coord = (y_atom >= 0) ? X[y_atom] : invalid;

    bool x_excluded = (x_atom >= 0) ? D_excluded[x_atom] : true;
    bool y_excluded = (y_atom >= 0) ? D_excluded[y_atom] : true;

    double x_charge = (x_atom >= 0) ? D_ccharges[D_charges[x_atom].code - 1].charge : 0.0;
    double y_charge = (y_atom >= 0) ? D_ccharges[D_charges[y_atom].code - 1].charge : 0.0;

    catype_t x_type = (x_atom >= 0) ? D_catypes[D_atypes[x_atom].code - 1] : catype_t{};
    catype_t y_type = (y_atom >= 0) ? D_catypes[D_atypes[y_atom].code - 1] : catype_t{};
    double x_aii_normal = x_type.aii_normal;
    double y_aii_normal = y_type.aii_normal;
    double x_aii_14 = x_type.aii_1_4;
    double y_aii_14 = y_type.aii_1_4;
    double x_bii_normal = x_type.bii_normal;
    double y_bii_normal = y_type.bii_normal;
    double x_bii_14 = x_type.bii_1_4;
    double y_bii_14 = y_type.bii_1_4;

    double3 x_force = {0.0, 0.0, 0.0};
    double3 y_force = {0.0, 0.0, 0.0};

    double evdw_sum = 0;
    double ecoul_sum = 0;

    const unsigned mask = 0xffffffffu;

    auto is_valid = [&]() -> bool {
        if (x_idx >= N || y_idx >= N) {
            return false;
        }
        if (tile_x == tile_y && x_idx <= y_idx) {
            return false;
        }
        if (x_excluded || y_excluded) {
            return false;
        }
        bool bond23 = (D_LJ_matrix[x_atom * n_atoms_solute + y_atom] == 3);
        if (bond23) {
            return false;
        }

        return true;
    };

    // shuffle to get next y
    auto do_shuffle = [&]() {
        const int src = (lane + 1) & 31;
        y_idx = __shfl_sync(mask, y_idx, src);
        y_atom = __shfl_sync(mask, y_atom, src);
        y_coord = shfl_coord(y_coord, src, mask);
        y_excluded = __shfl_sync(mask, y_excluded, src);
        y_charge = __shfl_sync(mask, y_charge, src);
        y_aii_normal = __shfl_sync(mask, y_aii_normal, src);
        y_aii_14 = __shfl_sync(mask, y_aii_14, src);
        y_bii_normal = __shfl_sync(mask, y_bii_normal, src);
        y_bii_14 = __shfl_sync(mask, y_bii_14, src);

        // shuffle y_force
        y_force.x = shfl(y_force.x, src, mask);
        y_force.y = shfl(y_force.y, src, mask);
        y_force.z = shfl(y_force.z, src, mask);
    };

    for (int i = 0; i < 32; i++) {
        if (is_valid()) {
            bool bond14 = (D_LJ_matrix[y_atom * n_atoms_solute + x_atom] == 1);
            double scaling = bond14 ? D_topo.el14_scale : 1.0;
            double ai_aii = bond14 ? x_aii_14 : x_aii_normal;
            double aj_aii = bond14 ? y_aii_14 : y_aii_normal;
            double ai_bii = bond14 ? x_bii_14 : x_bii_normal;
            double aj_bii = bond14 ? y_bii_14 : y_bii_normal;
            double evdw = 0.0, ecoul = 0.0, dv = 0.0;
            calculate_unforce_bound(
                x_coord, y_coord, x_charge, y_charge, ai_aii, aj_aii, ai_bii, aj_bii,
                D_topo.coulomb_constant, scaling, evdw, ecoul, dv);
            evdw_sum += evdw;
            ecoul_sum += ecoul;
            double3 d = {y_coord.x - x_coord.x, y_coord.y - x_coord.y, y_coord.z - x_coord.z};
            x_force.x -= dv * d.x;
            x_force.y -= dv * d.y;
            x_force.z -= dv * d.z;

            y_force.x += dv * d.x;
            y_force.y += dv * d.y;
            y_force.z += dv * d.z;
        }
        do_shuffle();
    }

    if (x_atom >= 0) {
        atomicAdd(&DV_X[x_atom].x, x_force.x);
        atomicAdd(&DV_X[x_atom].y, x_force.y);
        atomicAdd(&DV_X[x_atom].z, x_force.z);
    }
    if (y_atom >= 0) {
        atomicAdd(&DV_X[y_atom].x, y_force.x);
        atomicAdd(&DV_X[y_atom].y, y_force.y);
        atomicAdd(&DV_X[y_atom].z, y_force.z);
    }

    for (int offset = 16; offset > 0; offset >>= 1) {
        evdw_sum += __shfl_down_sync(0xffffffffu, evdw_sum, offset);
        ecoul_sum += __shfl_down_sync(0xffffffffu, ecoul_sum, offset);
    }
    if (lane == 0) {
        atomicAdd(Evdw_TOT, evdw_sum);
        atomicAdd(Ecoul_TOT, ecoul_sum);
    }
}

__global__ void calc_pp_dvel_matrix(int n_patoms, int n_atoms_solute,
                                    coord_t* X, double* Evdw, double* Ecoul, dvel_t* PP_MAT,
                                    ccharge_t* D_ccharges, charge_t* D_charges, catype_t* D_catypes, atype_t* D_atypes, p_atom_t* D_patoms, int* D_LJ_matrix, bool* D_excluded, topo_t D_topo) {
    // Block index
    int bx = blockIdx.x;
    int by = blockIdx.y;

    // Thread index
    int tx = threadIdx.x;
    int ty = threadIdx.y;

    int aStart = BLOCK_SIZE * by;
    int bStart = BLOCK_SIZE * bx;

    __shared__ double Ecoul_S[BLOCK_SIZE][BLOCK_SIZE];
    __shared__ double Evdw_S[BLOCK_SIZE][BLOCK_SIZE];

    Ecoul_S[ty][tx] = 0;
    Evdw_S[ty][tx] = 0;

    if (aStart + ty >= n_patoms) return;
    if (bStart + tx >= n_patoms) return;

    __shared__ coord_t Xs[BLOCK_SIZE];
    __shared__ coord_t Ys[BLOCK_SIZE];
    __shared__ int LJs[BLOCK_SIZE * BLOCK_SIZE];
    __shared__ bool excluded_s[2 * BLOCK_SIZE];

    int pi = D_patoms[aStart + ty].a - 1;
    int pj = D_patoms[bStart + tx].a - 1;

    if (tx == 0) {
        Xs[ty] = X[pi];
        excluded_s[ty] = D_excluded[pi];
    }

    if (ty == 0) {
        Ys[tx] = X[pj];
        excluded_s[BLOCK_SIZE + tx] = D_excluded[pj];
    }
    LJs[ty * BLOCK_SIZE + tx] = D_LJ_matrix[pi * n_atoms_solute + pj];

    __syncthreads();

    if (bx < by || (bx == by && tx < ty)) return;

    dvel_t patom_a, patom_b;
    memset(&patom_a, 0, sizeof(dvel_t));
    memset(&patom_b, 0, sizeof(dvel_t));

    int row = by * BLOCK_SIZE + ty;
    int column = bx * BLOCK_SIZE + tx;

    if (bx != by || tx != ty) {
        double evdw = 0, ecoul = 0;
        calc_pp_dvel_matrix_incr(ty, pi, tx, pj, Xs, Ys, LJs, excluded_s, &evdw, &ecoul, &patom_a,
                                 &patom_b, D_ccharges, D_charges, D_catypes, D_atypes, D_patoms, D_topo);
        Evdw_S[ty][tx] = evdw;
        Ecoul_S[ty][tx] = ecoul;
    }

    PP_MAT[row * n_patoms + column] = patom_a;
    PP_MAT[column * n_patoms + row] = patom_b;

    __syncthreads();

    int rowlen = (n_patoms + BLOCK_SIZE - 1) / BLOCK_SIZE;

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

        // printf("bx = %d by = %d tot_Evdw = %f tot_Ecoul = %f\n", bx, by, Evdw[rowlen * by + bx], Ecoul[rowlen * by + bx]);
    }

    __syncthreads();
}

__global__ void calc_pp_dvel_vector(int n_patoms, dvel_t* DV_X, dvel_t* PP_MAT, p_atom_t* D_patoms) {
    int row = blockIdx.x * blockDim.x + threadIdx.x;
    if (row >= n_patoms) return;

    dvel_t dP;
    dP.x = 0;
    dP.y = 0;
    dP.z = 0;

    for (int i = 0; i < n_patoms; i++) {
        if (i != row) {
            dP.x += PP_MAT[i + n_patoms * row].x;
            dP.y += PP_MAT[i + n_patoms * row].y;
            dP.z += PP_MAT[i + n_patoms * row].z;
        }
    }

    int p = D_patoms[row].a - 1;

    DV_X[p].x += dP.x;
    DV_X[p].y += dP.y;
    DV_X[p].z += dP.z;

    __syncthreads();
}

__global__ void calc_energy_sum(int rows, int columns, double* Evdw_TOT, double* Ecoul_TOT, double* Evdw, double* Ecoul, bool upper_diagonal) {
    int tx = threadIdx.x;
    int ty = threadIdx.y;

    __shared__ double Ecoul_S[BLOCK_SIZE][BLOCK_SIZE];
    __shared__ double Evdw_S[BLOCK_SIZE][BLOCK_SIZE];

    // TODO: better way to distribute upper diagonal over threads? Seems like threads in left bottom corner have less work
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
}  // namespace CudaNonbondedPPForce
void calc_nonbonded_pp_forces_host_v2() {
    using namespace CudaNonbondedPPForce;

    int mem_size_DV_X = n_atoms_solute * sizeof(dvel_t);

    CudaContext& ctx = CudaContext::instance();
    auto X = ctx.d_coords;
    auto DV_X = ctx.d_dvelocities;
    auto D_ccharges = ctx.d_ccharges;
    auto D_charges = ctx.d_charges;
    auto D_catypes = ctx.d_catypes;
    auto D_atypes = ctx.d_atypes;
    auto D_patoms = ctx.d_p_atoms;
    auto D_LJ_matrix = ctx.d_LJ_matrix;
    auto D_excluded = ctx.d_excluded;

    cudaMemset(D_PP_evdw_TOT, 0, sizeof(double));
    cudaMemset(D_PP_ecoul_TOT, 0, sizeof(double));

    const int thread_num = 128;
    dim3 block_sz = dim3(thread_num);

    int tile_num_per_block = thread_num / 32;
    int nBlocks = (n_patoms + 31) / 32;
    int total_tiles = nBlocks * (nBlocks + 1) / 2;

    int grid_sz = (total_tiles + tile_num_per_block - 1) / tile_num_per_block;
    dim3 grid = dim3(grid_sz);

    calc_pp<<<grid, block_sz>>>(n_patoms, n_atoms_solute, D_charges, D_ccharges,
                                D_atypes, D_catypes, topo, D_excluded,
                                D_LJ_matrix, D_patoms, X, DV_X, D_PP_evdw_TOT,
                                D_PP_ecoul_TOT);

    cudaDeviceSynchronize();
    cudaMemcpy(dvelocities, DV_X, mem_size_DV_X, cudaMemcpyDeviceToHost);

    double PP_evdw_TOT, PP_ecoul_TOT;
    cudaMemcpy(&PP_evdw_TOT, D_PP_evdw_TOT, sizeof(double), cudaMemcpyDeviceToHost);
    cudaMemcpy(&PP_ecoul_TOT, D_PP_ecoul_TOT, sizeof(double), cudaMemcpyDeviceToHost);

    E_nonbond_pp.Uvdw += PP_evdw_TOT;
    E_nonbond_pp.Ucoul += PP_ecoul_TOT;
}

void init_nonbonded_pp_force_kernel_data() {
    using namespace CudaNonbondedPPForce;

    if (!is_initialized) {
        check_cudaMalloc((void**)&D_PP_evdw_TOT, sizeof(double));
        check_cudaMalloc((void**)&D_PP_ecoul_TOT, sizeof(double));

        is_initialized = true;
    }
}

void cleanup_nonbonded_pp_force() {
    using namespace CudaNonbondedPPForce;
    if (is_initialized) {
        cudaFree(D_PP_evdw_TOT);
        cudaFree(D_PP_ecoul_TOT);
        is_initialized = false;
    }
}
