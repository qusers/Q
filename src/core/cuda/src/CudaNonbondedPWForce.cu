#include "cuda/include/CudaContext.cuh"
#include "cuda/include/CudaNonbondedPWForce.cuh"

namespace CudaNonbondedPWForce {
// Declare any necessary static variables or device pointers here
bool is_initialized = false;
struct calc_pw_t {
    dvel_t P;
    dvel_t W;
};

calc_pw_t *PW_MAT, *h_PW_MAT;
double *D_PW_Evdw, *D_PW_Ecoul;
double *D_PW_evdw_TOT, *D_PW_ecoul_TOT, PW_evdw_TOT, PW_ecoul_TOT;

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

__device__ void calc_pw_dvel_matrix_incr(int row, int pi, int column, int j, int n_atoms_solute,
                                         coord_t* Xs, coord_t* Ws, bool* excluded_s, double* Evdw, double* Ecoul, calc_pw_t* pw,
                                         ccharge_t* D_ccharges, charge_t* D_charges, catype_t* D_catypes, atype_t* D_atypes, p_atom_t* D_patoms, topo_t D_topo) {
    coord_t da;
    double r2a, ra, r6a;
    double Vela, V_a, V_b;
    double dva;
    double qi, qj;
    double ai_aii, aj_aii, ai_bii, aj_bii;
    catype_t ai_type, aj_type;
    atype_t i_type, j_type;

    if (excluded_s[row]) return;
    qi = D_ccharges[D_charges[pi].code - 1].charge;
    qj = D_ccharges[D_charges[n_atoms_solute + j].code - 1].charge;  // TODO: FIX THIS!!! WILL NOT WORK WITH QATOMS!!!!!

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
    //     printf("Vela = %f V_a = %f V_b = %f P = %f %f %f ai_aii = %f aj_aii = %f\n", Vela, V_a, V_b, pw->P.x, pw->P.y, pw->P.z, ai_aii, aj_aii);
    // }

    // if (pi < 100 && j < 100) printf("Evdw = %f Ecoul = %f\n", *Evdw, *Ecoul);
}

__global__ void calc_pw_dvel_vector_row(int n_patoms, int n_waters, dvel_t* DV_X, dvel_t* DV_W, calc_pw_t* PW_MAT, p_atom_t* D_patoms) {
    int row = blockIdx.x * blockDim.x + threadIdx.x;
    if (row >= n_patoms) return;

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

__global__ void calc_pw_dvel_vector_column(int n_patoms, int n_waters, dvel_t* DV_X, dvel_t* DV_W, calc_pw_t* PW_MAT) {
    int column = blockIdx.x * blockDim.x + threadIdx.x;
    if (column >= 3 * n_waters) return;

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

__global__ void calc_pw_dvel_matrix(int n_patoms, int n_atoms_solute, int n_waters,
                                    coord_t* X, coord_t* W, double* Evdw, double* Ecoul, calc_pw_t* PW_MAT,
                                    ccharge_t* D_ccharges, charge_t* D_charges, catype_t* D_catypes, atype_t* D_atypes, p_atom_t* D_patoms, bool* D_excluded, topo_t D_topo) {
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

    // if (bx == 0 && by == 0) printf("bx = %d by = %d tx = %d ty = %d\n", bx, by, tx, ty);

    int aStart = BLOCK_SIZE * by;
    int bStart = BLOCK_SIZE * bx;

    if (aStart + ty >= n_patoms) return;
    if (bStart + tx >= 3 * n_waters) return;

    // if (bx == 8 && by == 1) printf("bx = %d by = %d tx = %d ty = %d\n", bx, by, tx, ty);

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

    //     printf("Ys[%d] = %f Xs[%d] = %f\n", 3 * ty, Ys[3 * ty], 3 * tx, Xs[3 * tx]);
    // }

    // if (bx == 8 && by == 1) printf("bx = %d by = %d tx = %d ty = %d\n", bx, by, tx, ty);
    // __device__ void calc_pw_dvel_matrix_incr(int row, int pi, int column, int j, int n_patoms,
    // coord_t *Ps, coord_t *Xs, double *Evdw, double *Ecoul, calc_pw_t *pw,
    // ccharge_t *D_ccharges, charge_t *D_charges, catype_t *D_catypes, atype_t *D_atypes, p_atom_t *D_patoms)
    double evdw = 0, ecoul = 0;
    calc_pw_dvel_matrix_incr(ty, pi, tx, bStart + tx, n_atoms_solute, Xs, Ws, excluded_s, &evdw, &ecoul, &pw, D_ccharges, D_charges, D_catypes, D_atypes, D_patoms, D_topo);
    Evdw_S[ty][tx] = evdw;
    Ecoul_S[ty][tx] = ecoul;

    // if (row == 0 && column == 1) {
    //     printf("water_a = %f %f %f water_b = %f %f %f\n", water_a[0].x, water_a[0].y, water_a[0].z, water_b[0].x, water_b[0].y, water_b[0].z);
    // }

    // if (bx == 8 && by == 1) printf("n_qatoms = %d\n", n_qatoms);
    // if (bx == 8 && by == 1) printf("qi = %d j = %d charge[%d] = %f\n", row, column, row + n_qatoms, D_qcharges[row + n_qatoms * 1].q);

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

}  // namespace CudaNonbondedPWForce

void calc_nonbonded_pw_forces_host_v2() {
    using namespace CudaNonbondedPWForce;
    int mem_size_W = 3 * n_waters * sizeof(coord_t);
    int mem_size_DV_W = 3 * n_waters * sizeof(dvel_t);
    int mem_size_DV_X = n_atoms_solute * sizeof(dvel_t);
    int mem_size_PW_MAT = 3 * n_waters * n_patoms * sizeof(calc_pw_t);

    int n_blocks_p = (n_patoms + BLOCK_SIZE - 1) / BLOCK_SIZE;
    int n_blocks_w = (3 * n_waters + BLOCK_SIZE - 1) / BLOCK_SIZE;

    int mem_size_PW_Evdw = n_blocks_p * n_blocks_w * sizeof(double);
    int mem_size_PW_Ecoul = n_blocks_p * n_blocks_w * sizeof(double);

    if (!is_initialized) {
#ifdef DEBUG
        printf("Allocating PW_MAT\n");
#endif
        check_cudaMalloc((void**)&PW_MAT, mem_size_PW_MAT);

#ifdef DEBUG
        printf("Allocating D_PW_Evdw\n");
#endif
        check_cudaMalloc((void**)&D_PW_Evdw, mem_size_PW_Evdw);
#ifdef DEBUG
        printf("Allocating D_PW_Ecoul\n");
#endif
        check_cudaMalloc((void**)&D_PW_Ecoul, mem_size_PW_Ecoul);

        check_cudaMalloc((void**)&D_PW_evdw_TOT, sizeof(double));
        check_cudaMalloc((void**)&D_PW_ecoul_TOT, sizeof(double));

#ifdef DEBUG
        printf("All GPU solvent memory allocated\n");
#endif

        h_PW_MAT = (calc_pw_t*)malloc(mem_size_PW_MAT);
        is_initialized = true;
    }

    CudaContext& ctx = CudaContext::instance();
    auto X = ctx.d_coords;
    auto DV_X = ctx.d_dvelocities;
    auto D_ccharges = ctx.d_ccharges;
    auto D_charges = ctx.d_charges;
    auto D_catypes = ctx.d_catypes;
    auto D_atypes = ctx.d_atypes;
    auto D_patoms = ctx.d_p_atoms;
    auto D_excluded = ctx.d_excluded;

    dim3 threads, grid;

    threads = dim3(BLOCK_SIZE, BLOCK_SIZE);
    grid = dim3((3 * n_waters + BLOCK_SIZE - 1) / threads.x, (n_patoms + BLOCK_SIZE - 1) / threads.y);

    calc_pw_dvel_matrix<<<grid, threads>>>(n_patoms, n_atoms_solute, n_waters, X, X + n_atoms_solute, D_PW_Evdw, D_PW_Ecoul, PW_MAT,
                                           D_ccharges, D_charges, D_catypes, D_atypes, D_patoms, D_excluded, topo);
    calc_pw_dvel_vector_column<<<((3 * n_waters + BLOCK_SIZE - 1) / BLOCK_SIZE), BLOCK_SIZE>>>(n_patoms, n_waters, DV_X, DV_X + n_atoms_solute, PW_MAT);
    calc_pw_dvel_vector_row<<<((n_patoms + BLOCK_SIZE - 1) / BLOCK_SIZE), BLOCK_SIZE>>>(n_patoms, n_waters, DV_X, DV_X + n_atoms_solute, PW_MAT, D_patoms);

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
                printf("PW_MAT[%d][%d].P = %f %f %f\n", i, j, h_PW_MAT[3 * i * n_waters + j].P.x, h_PW_MAT[3 * i * n_waters + j].P.y, h_PW_MAT[3 * i * n_waters + j].P.z);
        }
    }
#endif

    cudaMemcpy(dvelocities, DV_X, mem_size_DV_X, cudaMemcpyDeviceToHost);
    cudaMemcpy(&dvelocities[n_atoms_solute], DV_X + n_atoms_solute, mem_size_DV_W, cudaMemcpyDeviceToHost);

    calc_energy_sum<<<1, threads>>>(n_blocks_p, n_blocks_w, D_PW_evdw_TOT, D_PW_ecoul_TOT, D_PW_Evdw, D_PW_Ecoul, false);

    cudaMemcpy(&PW_evdw_TOT, D_PW_evdw_TOT, sizeof(double), cudaMemcpyDeviceToHost);
    cudaMemcpy(&PW_ecoul_TOT, D_PW_ecoul_TOT, sizeof(double), cudaMemcpyDeviceToHost);

    E_nonbond_pw.Uvdw += PW_evdw_TOT;
    E_nonbond_pw.Ucoul += PW_ecoul_TOT;

    // for (int i = 0; i < n_atoms; i++) {
    //     printf("dvelocities[%d] = %f %f %f\n", i, dvelocities[i].x, dvelocities[i].y, dvelocities[i].z);
    // }
}

void cleanup_nonbonded_pw_force() {
    using namespace CudaNonbondedPWForce;
    if (is_initialized) {
        cudaFree(PW_MAT);
        cudaFree(D_PW_Evdw);
        cudaFree(D_PW_Ecoul);
        cudaFree(D_PW_evdw_TOT);
        cudaFree(D_PW_ecoul_TOT);
        free(h_PW_MAT);
        is_initialized = false;
    }
}