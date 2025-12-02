#include <iostream>

#include "cuda/include/CudaContext.cuh"
#include "cuda/include/CudaNonbondedQPForce.cuh"
#include "system.h"

namespace CudaNonbondedQPForce {
struct calc_qw_t {
    dvel_t Q;
    dvel_t O;
    dvel_t H1;
    dvel_t H2;
};

struct calc_qp_t {
    dvel_t Q;
    dvel_t P;
};
bool is_initialized = false;
double *D_QP_Evdw, *D_QP_Ecoul, *h_QP_Evdw, *h_QP_Ecoul;
calc_qp_t *QP_MAT, *h_QP_MAT;
double *D_QP_evdw_TOT, *D_QP_ecoul_TOT, QP_evdw_TOT, QP_ecoul_TOT;

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
__global__ void calc_qp_dvel_vector_row(int n_qatoms, int n_patoms, dvel_t* DV_X, calc_qp_t* QP_MAT, q_atom_t* D_qatoms) {
    int row = blockIdx.x * blockDim.x + threadIdx.x;
    if (row >= n_qatoms) return;

    dvel_t dQ;

    dQ.x = 0;
    dQ.y = 0;
    dQ.z = 0;

    for (int i = 0; i < n_patoms; i++) {
        dQ.x += QP_MAT[i + n_patoms * row].Q.x;
        dQ.y += QP_MAT[i + n_patoms * row].Q.y;
        dQ.z += QP_MAT[i + n_patoms * row].Q.z;
    }

    int q = D_qatoms[row].a - 1;

    DV_X[q].x += dQ.x;
    DV_X[q].y += dQ.y;
    DV_X[q].z += dQ.z;

    __syncthreads();
}

__global__ void calc_qp_dvel_vector_column(int n_qatoms, int n_patoms, dvel_t* DV_X, calc_qp_t* QP_MAT, p_atom_t* D_patoms) {
    int column = blockIdx.x * blockDim.x + threadIdx.x;
    if (column >= n_patoms) return;

    dvel_t dP;

    dP.x = 0;
    dP.y = 0;
    dP.z = 0;

    for (int i = 0; i < n_qatoms; i++) {
        dP.x += QP_MAT[column + n_patoms * i].P.x;
        dP.y += QP_MAT[column + n_patoms * i].P.y;
        dP.z += QP_MAT[column + n_patoms * i].P.z;
    }

    int p = D_patoms[column].a - 1;

    DV_X[p].x += dP.x;
    DV_X[p].y += dP.y;
    DV_X[p].z += dP.z;

    __syncthreads();
}

__device__ void calc_qp_dvel_matrix_incr(int row, int qi, int column, int pj, int n_lambdas, int n_qatoms,
                                         coord_t* Qs, coord_t* Ps, int* LJs, bool* excluded_s, double Evdw_S[BLOCK_SIZE][2 * BLOCK_SIZE], double Ecoul_S[BLOCK_SIZE][2 * BLOCK_SIZE], calc_qp_t* qp,
                                         q_catype_t* D_qcatypes, q_atype_t* D_qatypes, q_charge_t* D_qcharges, p_atom_t* D_patoms, q_atom_t* D_qatoms, double* D_lambdas,
                                         catype_t* D_catypes, atype_t* D_atypes, ccharge_t* D_ccharges, charge_t* D_charges, topo_t D_topo) {
    coord_t da;
    double r2, r6, r;
    double ai_aii, aj_aii, ai_bii, aj_bii;
    q_catype_t qi_type;
    catype_t aj_type;
    bool bond23, bond14;
    double scaling, Vel, V_a, V_b, dv;

    int j = D_patoms[pj].a - 1;

    bond23 = LJs[row * BLOCK_SIZE + column] == 3;
    bond14 = LJs[row * BLOCK_SIZE + column] == 1;

    if (bond23) return;
    if (excluded_s[row] || excluded_s[BLOCK_SIZE + column]) return;

    scaling = bond14 ? D_topo.el14_scale : 1;

    da.x = Qs[row].x - Ps[column].x;
    da.y = Qs[row].y - Ps[column].y;
    da.z = Qs[row].z - Ps[column].z;

    r2 = pow(da.x, 2) + pow(da.y, 2) + pow(da.z, 2);

    r6 = r2 * r2 * r2;
    r2 = 1 / r2;
    r = sqrt(r2);

    for (int state = 0; state < n_lambdas; state++) {
        qi_type = D_qcatypes[D_qatypes[qi + n_qatoms * state].code - 1];
        aj_type = D_catypes[D_atypes[j].code - 1];

        ai_aii = bond14 ? qi_type.Ai_14 : qi_type.Ai;
        aj_aii = bond14 ? aj_type.aii_1_4 : aj_type.aii_normal;
        ai_bii = bond14 ? qi_type.Bi_14 : qi_type.Bi;
        aj_bii = bond14 ? aj_type.bii_1_4 : aj_type.bii_normal;

        Vel = D_topo.coulomb_constant * scaling * D_qcharges[qi + n_qatoms * state].q * D_ccharges[D_charges[j].code - 1].charge * r;
        V_a = ai_aii * aj_aii / (r6 * r6);
        V_b = ai_bii * aj_bii / r6;
        dv = r2 * (-Vel - (12 * V_a - 6 * V_b)) * D_lambdas[state];

        // if (state == 0 && qi == 0 && pj == 1) {
        //     printf("crg_q = %f crg_j = %f r = %f\n", D_qcharges[qi + n_qatoms * state].q, D_ccharges[D_charges[pj].code - 1].charge, r);
        //     printf("ai_aii = %f aj_aii = %f ai_bii = %f aj_bii = %f\n", ai_aii, aj_aii, ai_bii, aj_bii);
        // }

        // Update forces
        qp->Q.x += dv * da.x;
        qp->Q.y += dv * da.y;
        qp->Q.z += dv * da.z;
        qp->P.x -= dv * da.x;
        qp->P.y -= dv * da.y;
        qp->P.z -= dv * da.z;

        // Update Q totals
        Ecoul_S[row][state * BLOCK_SIZE + column] += Vel;
        Evdw_S[row][state * BLOCK_SIZE + column] += (V_a - V_b);
    }
}

__global__ void calc_qp_dvel_matrix(int n_qatoms, int n_patoms, int n_lambdas, int n_atoms_solute,
                                    coord_t* X, double* Evdw, double* Ecoul, calc_qp_t* QP_MAT,
                                    q_catype_t* D_qcatypes, q_atype_t* D_qatypes, q_charge_t* D_qcharges, p_atom_t* D_patoms, q_atom_t* D_qatoms, double* D_lambdas,
                                    int* D_LJ_matrix, bool* D_excluded, catype_t* D_catypes, atype_t* D_atypes, ccharge_t* D_ccharges, charge_t* D_charges, topo_t D_topo) {
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
    int bStart = BLOCK_SIZE * bx;

    if (aStart + ty >= n_qatoms) return;
    if (bStart + tx >= n_patoms) return;

    int qi = D_qatoms[aStart + ty].a - 1;
    int pj = D_patoms[bStart + tx].a - 1;

    __shared__ coord_t Qs[BLOCK_SIZE];
    __shared__ coord_t Ps[BLOCK_SIZE];
    __shared__ int LJs[BLOCK_SIZE * BLOCK_SIZE];
    __shared__ bool excluded_s[2 * BLOCK_SIZE];

    if (tx == 0) {
        Qs[ty] = X[qi];
        excluded_s[ty] = D_excluded[qi];
    }

    if (ty == 0) {
        Ps[tx] = X[pj];
        excluded_s[BLOCK_SIZE + tx] = D_excluded[pj];
    }
    LJs[ty * BLOCK_SIZE + tx] = D_LJ_matrix[qi * n_atoms_solute + pj];

    __syncthreads();

    calc_qp_t qp;
    memset(&qp, 0, sizeof(calc_qw_t));

    int row = by * BLOCK_SIZE + ty;
    int column = bx * BLOCK_SIZE + tx;

    calc_qp_dvel_matrix_incr(ty, row, tx, column, n_lambdas, n_qatoms, Qs, Ps, LJs, excluded_s, Evdw_S, Ecoul_S,
                             &qp, D_qcatypes, D_qatypes, D_qcharges, D_patoms, D_qatoms, D_lambdas, D_catypes, D_atypes, D_ccharges, D_charges, D_topo);

    QP_MAT[n_patoms * row + column] = qp;

    __syncthreads();

    int rowlen = (n_patoms + BLOCK_SIZE - 1) / BLOCK_SIZE;
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

}  // namespace CudaNonbondedQPForce

void calc_nonbonded_qp_forces_host_v2() {
    using namespace CudaNonbondedQPForce;

    int n_blocks_q = (n_qatoms + BLOCK_SIZE - 1) / BLOCK_SIZE;
    int n_blocks_p = (n_patoms + BLOCK_SIZE - 1) / BLOCK_SIZE;

    CudaContext& ctx = CudaContext::instance();
    auto X = ctx.d_coords;
    auto DV_X = ctx.d_dvelocities;
    auto D_qcatypes = ctx.d_q_catypes;
    auto D_qatypes = ctx.d_q_atypes;
    auto D_qcharges = ctx.d_q_charges;
    auto D_patoms = ctx.d_p_atoms;
    auto D_qatoms = ctx.d_q_atoms;
    auto D_lambdas = ctx.d_lambdas;
    auto D_LJ_matrix = ctx.d_LJ_matrix;
    auto D_excluded = ctx.d_excluded;
    auto D_catypes = ctx.d_catypes;
    auto D_atypes = ctx.d_atypes;
    auto D_ccharges = ctx.d_ccharges;
    auto D_charges = ctx.d_charges;

    dim3 threads, grid;

    threads = dim3(BLOCK_SIZE, BLOCK_SIZE);
    grid = dim3((n_patoms + BLOCK_SIZE - 1) / threads.x, (n_qatoms + BLOCK_SIZE - 1) / threads.y);

    calc_qp_dvel_matrix<<<grid, threads>>>(n_qatoms, n_patoms, n_lambdas, n_atoms_solute, X, D_QP_Evdw, D_QP_Ecoul, QP_MAT,
                                           D_qcatypes, D_qatypes, D_qcharges, D_patoms, D_qatoms, D_lambdas, D_LJ_matrix, D_excluded,
                                           D_catypes, D_atypes, D_ccharges, D_charges, topo);
    calc_qp_dvel_vector_column<<<((n_patoms + BLOCK_SIZE - 1) / BLOCK_SIZE), BLOCK_SIZE>>>(n_qatoms, n_patoms, DV_X, QP_MAT, D_patoms);
    calc_qp_dvel_vector_row<<<((n_qatoms + BLOCK_SIZE - 1) / BLOCK_SIZE), BLOCK_SIZE>>>(n_qatoms, n_patoms, DV_X, QP_MAT, D_qatoms);

#ifdef DEBUG
    cudaMemcpy(h_QP_MAT, QP_MAT, mem_size_QP_MAT, cudaMemcpyDeviceToHost);
#endif

#ifdef DEBUG
    for (int i = 0; i < n_qatoms; i++) {
        for (int j = 0; j < n_patoms; j++) {
            if (i == 0)
                // if (h_QP_MAT[i * n_patoms + j].Q.x > 100)
                printf("QP_MAT[%d][%d].Q = %f %f %f\n", i, j, h_QP_MAT[i * n_patoms + j].Q.x, h_QP_MAT[i * n_patoms + j].Q.y, h_QP_MAT[i * n_patoms + j].Q.z);
        }
    }
#endif
    int mem_size_DV_X = n_atoms * sizeof(dvel_t);
    cudaMemcpy(dvelocities, DV_X, mem_size_DV_X, cudaMemcpyDeviceToHost);

    // TODO: make Evdw & Ecoul work for # of states > 2
    for (int state = 0; state < min(2, n_lambdas); state++) {
        calc_energy_sum<<<1, threads>>>(n_blocks_q, n_blocks_p, D_QP_evdw_TOT, D_QP_ecoul_TOT, &D_QP_Evdw[state * n_blocks_p * n_blocks_q], &D_QP_Ecoul[state * n_blocks_p * n_blocks_q], false);

        cudaMemcpy(&QP_evdw_TOT, D_QP_evdw_TOT, sizeof(double), cudaMemcpyDeviceToHost);
        cudaMemcpy(&QP_ecoul_TOT, D_QP_ecoul_TOT, sizeof(double), cudaMemcpyDeviceToHost);

        EQ_nonbond_qp[state].Uvdw += QP_evdw_TOT;
        EQ_nonbond_qp[state].Ucoul += QP_ecoul_TOT;
    }
}

void init_nonbonded_qp_force_kernel_data() {
    using namespace CudaNonbondedQPForce;

    if (!is_initialized) {
        int n_blocks_q = (n_qatoms + BLOCK_SIZE - 1) / BLOCK_SIZE;
        int n_blocks_p = (n_patoms + BLOCK_SIZE - 1) / BLOCK_SIZE;
        // TODO make Evdw & Ecoul work for # of states > 2
        int mem_size_QP_Evdw = min(n_lambdas, 2) * n_blocks_q * n_blocks_p * sizeof(double);
        int mem_size_QP_Ecoul = min(n_lambdas, 2) * n_blocks_q * n_blocks_p * sizeof(double);
        int mem_size_QP_MAT = n_qatoms * n_patoms * sizeof(calc_qp_t);

        check_cudaMalloc((void**)&D_QP_Evdw, mem_size_QP_Evdw);
        check_cudaMalloc((void**)&D_QP_Ecoul, mem_size_QP_Ecoul);
        h_QP_Evdw = (double*)malloc(mem_size_QP_Evdw);
        h_QP_Ecoul = (double*)malloc(mem_size_QP_Ecoul);

        check_cudaMalloc((void**)&QP_MAT, mem_size_QP_MAT);
        h_QP_MAT = (calc_qp_t*)malloc(mem_size_QP_MAT);

        check_cudaMalloc((void**)&D_QP_evdw_TOT, sizeof(double));
        check_cudaMalloc((void**)&D_QP_ecoul_TOT, sizeof(double));

        is_initialized = true;
    }
}

void cleanup_nonbonded_qp_force() {
    using namespace CudaNonbondedQPForce;

    if (is_initialized) {
        cudaFree(D_QP_Evdw);
        cudaFree(D_QP_Ecoul);
        free(h_QP_Evdw);
        free(h_QP_Ecoul);
        cudaFree(QP_MAT);
        free(h_QP_MAT);
        cudaFree(D_QP_evdw_TOT);
        cudaFree(D_QP_ecoul_TOT);

        is_initialized = false;
    }
}