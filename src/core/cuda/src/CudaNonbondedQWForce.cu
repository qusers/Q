#include <iostream>

#include "cuda/include/CudaContext.cuh"
#include "cuda/include/CudaNonbondedQWForce.cuh"
#include "utils.h"

namespace CudaNonbondedQWForce {
struct calc_qw_t {
    dvel_t Q;
    dvel_t O;
    dvel_t H1;
    dvel_t H2;
};
bool is_initialized = false;
calc_qw_t *QW_MAT, *h_QW_MAT;

double *D_QW_Evdw, *D_QW_Ecoul, *h_QW_Evdw, *h_QW_Ecoul;
double *D_QW_evdw_TOT, *D_QW_ecoul_TOT, QW_evdw_TOT, QW_ecoul_TOT;

// General
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

__global__ void calc_qw_dvel_vector_row(int n_qatoms, int n_waters, dvel_t* DV_X, dvel_t* DV_W, calc_qw_t* MAT, q_atom_t* D_qatoms) {
    int row = blockIdx.x * blockDim.x + threadIdx.x;
    if (row >= n_qatoms) return;

    dvel_t dQ;

    dQ.x = 0;
    dQ.y = 0;
    dQ.z = 0;

    for (int i = 0; i < n_waters; i++) {
        dQ.x += MAT[i + n_waters * row].Q.x;
        dQ.y += MAT[i + n_waters * row].Q.y;
        dQ.z += MAT[i + n_waters * row].Q.z;
    }

    int q = D_qatoms[row].a - 1;

    DV_X[q].x += dQ.x;
    DV_X[q].y += dQ.y;
    DV_X[q].z += dQ.z;

    __syncthreads();
}

__global__ void calc_qw_dvel_vector_column(int n_qatoms, int n_waters, dvel_t* DV_X, dvel_t* DV_W, calc_qw_t* MAT) {
    int column = blockIdx.x * blockDim.x + threadIdx.x;
    if (column >= n_waters) return;

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

    for (int i = 0; i < n_qatoms; i++) {
        dO.x += MAT[column + n_waters * i].O.x;
        dO.y += MAT[column + n_waters * i].O.y;
        dO.z += MAT[column + n_waters * i].O.z;
        dH1.x += MAT[column + n_waters * i].H1.x;
        dH1.y += MAT[column + n_waters * i].H1.y;
        dH1.z += MAT[column + n_waters * i].H1.z;
        dH2.x += MAT[column + n_waters * i].H2.x;
        dH2.y += MAT[column + n_waters * i].H2.y;
        dH2.z += MAT[column + n_waters * i].H2.z;
    }

    DV_W[3 * column].x += dO.x;
    DV_W[3 * column].y += dO.y;
    DV_W[3 * column].z += dO.z;
    DV_W[3 * column + 1].x += dH1.x;
    DV_W[3 * column + 1].y += dH1.y;
    DV_W[3 * column + 1].z += dH1.z;
    DV_W[3 * column + 2].x += dH2.x;
    DV_W[3 * column + 2].y += dH2.y;
    DV_W[3 * column + 2].z += dH2.z;

    __syncthreads();
}

__device__ void calc_qw_dvel_matrix_incr(int row, int qi, int column, int n_lambdas, int n_qatoms, double crg_ow, double crg_hw, double A_O, double B_O,
                                         coord_t* Qs, coord_t* Ws, double Evdw_S[BLOCK_SIZE][2 * BLOCK_SIZE], double Ecoul_S[BLOCK_SIZE][2 * BLOCK_SIZE], calc_qw_t* qw,
                                         q_catype_t* D_qcatypes, q_atype_t* D_qatypes, q_charge_t* D_qcharges, q_atom_t* D_qatoms, double* D_lambdas, topo_t D_topo) {
    int j;
    coord_t dO, dH1, dH2;
    double r2O, rH1, rH2, r6O, rO, r2H1, r2H2;
    double dvO, dvH1, dvH2;
    double V_a, V_b, VelO, VelH1, VelH2;
    q_atype_t qa_type;
    q_catype_t qi_type;
    double ai_aii, ai_bii;

    j = 3 * column;
    dO.x = Ws[j].x - Qs[row].x;
    dO.y = Ws[j].y - Qs[row].y;
    dO.z = Ws[j].z - Qs[row].z;
    dH1.x = Ws[j + 1].x - Qs[row].x;
    dH1.y = Ws[j + 1].y - Qs[row].y;
    dH1.z = Ws[j + 1].z - Qs[row].z;
    dH2.x = Ws[j + 2].x - Qs[row].x;
    dH2.y = Ws[j + 2].y - Qs[row].y;
    dH2.z = Ws[j + 2].z - Qs[row].z;

    r2O = pow(dO.x, 2) + pow(dO.y, 2) + pow(dO.z, 2);
    rH1 = sqrt(1.0 / (pow(dH1.x, 2) + pow(dH1.y, 2) + pow(dH1.z, 2)));
    rH2 = sqrt(1.0 / (pow(dH2.x, 2) + pow(dH2.y, 2) + pow(dH2.z, 2)));
    r6O = r2O * r2O * r2O;
    r2O = 1.0 / r2O;
    rO = sqrt(r2O);
    r2H1 = rH1 * rH1;
    r2H2 = rH2 * rH2;

    // Reset potential
    dvO = 0;
    dvH1 = 0;
    dvH2 = 0;

    for (int state = 0; state < n_lambdas; state++) {
        qa_type = D_qatypes[qi + n_qatoms * state];
        qi_type = D_qcatypes[qa_type.code - 1];

        ai_aii = qi_type.Ai;
        ai_bii = qi_type.Bi;

        V_a = ai_aii * A_O / (r6O * r6O);
        V_b = ai_bii * B_O / (r6O);

        VelO = D_topo.coulomb_constant * crg_ow * D_qcharges[qi + n_qatoms * state].q * rO;
        VelH1 = D_topo.coulomb_constant * crg_hw * D_qcharges[qi + n_qatoms * state].q * rH1;
        VelH2 = D_topo.coulomb_constant * crg_hw * D_qcharges[qi + n_qatoms * state].q * rH2;

        dvO += r2O * (-VelO - (12 * V_a - 6 * V_b)) * D_lambdas[state];
        dvH1 -= r2H1 * VelH1 * D_lambdas[state];
        dvH2 -= r2H2 * VelH2 * D_lambdas[state];

        // Update Q totals
        Ecoul_S[row][state * BLOCK_SIZE + column] += (VelO + VelH1 + VelH2);
        Evdw_S[row][state * BLOCK_SIZE + column] += (V_a - V_b);
    }

    // Note r6O is not the usual 1/rO^6, but rather rO^6. be careful!!!

    // Update forces on Q-atom
    (*qw).Q.x -= (dvO * dO.x + dvH1 * dH1.x + dvH2 * dH2.x);
    (*qw).Q.y -= (dvO * dO.y + dvH1 * dH1.y + dvH2 * dH2.y);
    (*qw).Q.z -= (dvO * dO.z + dvH1 * dH1.z + dvH2 * dH2.z);

    // Update forces on water
    (*qw).O.x += dvO * dO.x;
    (*qw).O.y += dvO * dO.y;
    (*qw).O.z += dvO * dO.z;
    (*qw).H1.x += dvH1 * dH1.x;
    (*qw).H1.y += dvH1 * dH1.y;
    (*qw).H1.z += dvH1 * dH1.z;
    (*qw).H2.x += dvH2 * dH2.x;
    (*qw).H2.y += dvH2 * dH2.y;
    (*qw).H2.z += dvH2 * dH2.z;
}

__global__ void calc_qw_dvel_matrix(int n_qatoms, int n_waters, int n_lambdas, double crg_ow, double crg_hw, double A_O, double B_O,
                                    coord_t* X, coord_t* W, double* Evdw, double* Ecoul, calc_qw_t* MAT,
                                    q_catype_t* D_qcatypes, q_atype_t* D_qatypes, q_charge_t* D_qcharges, q_atom_t* D_qatoms, double* D_lambdas, topo_t D_topo) {
    // Block index
    int bx = blockIdx.x;
    int by = blockIdx.y;

    // Thread index
    int tx = threadIdx.x;
    int ty = threadIdx.y;

    // TODO implement >2 states on GPU
    __shared__ double Evdw_S[BLOCK_SIZE][2 * BLOCK_SIZE];
    __shared__ double Ecoul_S[BLOCK_SIZE][2 * BLOCK_SIZE];

    Ecoul_S[ty][tx] = 0;
    Evdw_S[ty][tx] = 0;
    Ecoul_S[ty][tx + BLOCK_SIZE] = 0;
    Evdw_S[ty][tx + BLOCK_SIZE] = 0;

    int aStart = BLOCK_SIZE * by;
    int bStart = 3 * BLOCK_SIZE * bx;

    if (aStart + ty >= n_qatoms) return;
    if (bStart + 3 * tx >= 3 * n_waters) return;

    __shared__ coord_t Qs[BLOCK_SIZE];
    __shared__ coord_t Ws[3 * BLOCK_SIZE];

    if (tx == 0) {
        Qs[ty] = X[D_qatoms[aStart + ty].a - 1];
    }

    if (ty == 0) {
        Ws[3 * tx] = W[bStart + 3 * tx];
        Ws[3 * tx + 1] = W[bStart + 3 * tx + 1];
        Ws[3 * tx + 2] = W[bStart + 3 * tx + 2];
    }

    __syncthreads();

    calc_qw_t qw;
    memset(&qw, 0, sizeof(calc_qw_t));

    int row = by * BLOCK_SIZE + ty;
    int column = bx * BLOCK_SIZE + tx;

    calc_qw_dvel_matrix_incr(ty, aStart + ty, tx, n_lambdas, n_qatoms, crg_ow, crg_hw, A_O, B_O, Qs, Ws, Evdw_S, Ecoul_S,
                             &qw, D_qcatypes, D_qatypes, D_qcharges, D_qatoms, D_lambdas, D_topo);

    MAT[column + n_waters * row] = qw;

    __syncthreads();

    int rowlen = (n_waters + BLOCK_SIZE - 1) / BLOCK_SIZE;
    int collen = (n_qatoms + BLOCK_SIZE - 1) / BLOCK_SIZE;

    if (tx == 0 && ty == 0) {
        // TODO implement >2 states on GPU
        for (int state = 0; state < min(2, n_lambdas); state++) {
            double tot_Evdw = 0;
            double tot_Ecoul = 0;
            for (int i = 0; i < BLOCK_SIZE; i++) {
                for (int j = 0; j < BLOCK_SIZE; j++) {
                    tot_Evdw += Evdw_S[i][j + state * BLOCK_SIZE];
                    tot_Ecoul += Ecoul_S[i][j + state * BLOCK_SIZE];
                }
            }
            Evdw[rowlen * collen * state + rowlen * by + bx] = tot_Evdw;
            Ecoul[rowlen * collen * state + rowlen * by + bx] = tot_Ecoul;
        }
    }

    __syncthreads();
}
}  // namespace CudaNonbondedQWForce

void calc_nonbonded_qw_forces_host_v2() {
    using namespace CudaNonbondedQWForce;

    int mem_size_X = n_atoms_solute * sizeof(coord_t);
    int mem_size_W = 3 * n_waters * sizeof(coord_t);
    int mem_size_DV_X = n_atoms_solute * sizeof(dvel_t);
    int mem_size_DV_W = 3 * n_waters * sizeof(dvel_t);
    int mem_size_MAT = 3 * n_waters * n_qatoms * sizeof(calc_qw_t);

    int n_blocks_q = (n_qatoms + BLOCK_SIZE - 1) / BLOCK_SIZE;
    int n_blocks_w = (n_waters + BLOCK_SIZE - 1) / BLOCK_SIZE;
    // TODO make Evdw & Ecoul work for # of states > 2
    int mem_size_QW_Evdw = min(n_lambdas, 2) * n_blocks_q * n_blocks_w * sizeof(double);
    int mem_size_QW_Ecoul = min(n_lambdas, 2) * n_blocks_q * n_blocks_w * sizeof(double);

    CudaContext& ctx = CudaContext::instance();
    auto X = ctx.d_coords;
    auto DV_X = ctx.d_dvelocities;
    auto D_qcatypes = ctx.d_q_catypes;
    auto D_qatypes = ctx.d_q_atypes;
    auto D_qcharges = ctx.d_q_charges;
    auto D_qatoms = ctx.d_q_atoms;
    auto D_lambdas = ctx.d_lambdas;

    dim3 threads, grid;

    threads = dim3(BLOCK_SIZE, BLOCK_SIZE);
    grid = dim3((n_waters + BLOCK_SIZE - 1) / threads.x, (n_qatoms + BLOCK_SIZE - 1) / threads.y);

    double evdw, ecoul;

    calc_qw_dvel_matrix<<<grid, threads>>>(n_qatoms, n_waters, n_lambdas, crg_ow, crg_hw, A_O, B_O, X, X + n_atoms_solute, D_QW_Evdw, D_QW_Ecoul,
                                           QW_MAT, D_qcatypes, D_qatypes, D_qcharges, D_qatoms, D_lambdas, topo);
    calc_qw_dvel_vector_column<<<((n_waters + BLOCK_SIZE - 1) / BLOCK_SIZE), BLOCK_SIZE>>>(n_qatoms, n_waters, DV_X, DV_X + n_atoms_solute, QW_MAT);
    calc_qw_dvel_vector_row<<<((n_qatoms + BLOCK_SIZE - 1) / BLOCK_SIZE), BLOCK_SIZE>>>(n_qatoms, n_waters, DV_X, DV_X + n_atoms_solute, QW_MAT, D_qatoms);

#ifdef DEBUG
    cudaMemcpy(h_QW_MAT, QW_MAT, mem_size_MAT, cudaMemcpyDeviceToHost);
#endif

#ifdef DEBUG
    for (int i = 0; i < n_qatoms; i++) {
        for (int j = 0; j < 3 * n_waters; j++) {
            if (h_QW_MAT[3 * i * n_waters + j].Q.x > 100)
                printf("QW_MAT[%d][%d].Q = %f %f %f\n", i, j, h_QW_MAT[i * 3 * n_waters + j].Q.x, h_QW_MAT[i * 3 * n_waters + j].Q.y, h_QW_MAT[i * 3 * n_waters + j].Q.z);
        }
    }
#endif

    cudaMemcpy(dvelocities, DV_X, mem_size_DV_X, cudaMemcpyDeviceToHost);
    cudaMemcpy(&dvelocities[n_atoms_solute], DV_X + n_atoms_solute, mem_size_DV_W, cudaMemcpyDeviceToHost);

    // TODO make Evdw & Ecoul work for # of states > 2
    for (int state = 0; state < min(2, n_lambdas); state++) {
        calc_energy_sum<<<1, threads>>>(n_blocks_q, n_blocks_w, D_QW_evdw_TOT, D_QW_ecoul_TOT, &D_QW_Evdw[state * n_blocks_w * n_blocks_q],
                                        &D_QW_Ecoul[state * n_blocks_w * n_blocks_q], false);

        cudaMemcpy(&QW_evdw_TOT, D_QW_evdw_TOT, sizeof(double), cudaMemcpyDeviceToHost);
        cudaMemcpy(&QW_ecoul_TOT, D_QW_ecoul_TOT, sizeof(double), cudaMemcpyDeviceToHost);

        EQ_nonbond_qw[state].Uvdw += QW_evdw_TOT;
        EQ_nonbond_qw[state].Ucoul += QW_ecoul_TOT;
    }
}

void init_nonbonded_qw_force_kernel_data() {
    using namespace CudaNonbondedQWForce;
    if (!is_initialized) {
        catype_t catype_ow;  // Atom type of first O atom

        catype_ow = catypes[atypes[n_atoms_solute].code - 1];

        A_O = catype_ow.aii_normal;
        B_O = catype_ow.bii_normal;

        int n_blocks_q = (n_qatoms + BLOCK_SIZE - 1) / BLOCK_SIZE;
        int n_blocks_w = (n_waters + BLOCK_SIZE - 1) / BLOCK_SIZE;
        int mem_size_MAT = 3 * n_waters * n_qatoms * sizeof(calc_qw_t);

        int mem_size_QW_Evdw = min(n_lambdas, 2) * n_blocks_q * n_blocks_w * sizeof(double);
        int mem_size_QW_Ecoul = min(n_lambdas, 2) * n_blocks_q * n_blocks_w * sizeof(double);
#ifdef DEBUG
        printf("Allocating QW_MAT\n");
#endif
        check_cudaMalloc((void**)&QW_MAT, mem_size_MAT);

#ifdef DEBUG
        printf("Allocating D_QW_Evdw\n");
#endif
        check_cudaMalloc((void**)&D_QW_Evdw, mem_size_QW_Evdw);
#ifdef DEBUG
        printf("Allocating D_QW_Ecoul\n");
#endif
        check_cudaMalloc((void**)&D_QW_Ecoul, mem_size_QW_Ecoul);

        check_cudaMalloc((void**)&D_QW_evdw_TOT, sizeof(double));
        check_cudaMalloc((void**)&D_QW_ecoul_TOT, sizeof(double));

        h_QW_Evdw = (double*)malloc(mem_size_QW_Evdw);
        h_QW_Ecoul = (double*)malloc(mem_size_QW_Ecoul);

        h_QW_MAT = (calc_qw_t*)malloc(mem_size_MAT);
        is_initialized = true;
    }
}

void cleanup_nonbonded_qw_force() {
    using namespace CudaNonbondedQWForce;

    if (is_initialized) {
        cudaFree(QW_MAT);
        cudaFree(D_QW_Evdw);
        cudaFree(D_QW_Ecoul);
        cudaFree(D_QW_evdw_TOT);
        cudaFree(D_QW_ecoul_TOT);

        free(h_QW_Evdw);
        free(h_QW_Ecoul);
        free(h_QW_MAT);

        is_initialized = false;
    }
}