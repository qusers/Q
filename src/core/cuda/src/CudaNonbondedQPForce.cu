#include <iostream>

#include "cuda/include/CudaContext.cuh"
#include "cuda/include/CudaNonbondedQPForce.cuh"
#include "system.h"

namespace CudaNonbondedQPForce {
bool is_initialized = false;
double *D_QP_Evdw, *D_QP_Ecoul, *h_QP_Evdw, *h_QP_Ecoul;
calc_qp_t *QP_MAT, *h_QP_MAT;

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

}  // namespace CudaNonbondedQPForce

void calc_nonbonded_qp_forces_host_v2() {
    using namespace CudaNonbondedQPForce;

    int n_blocks_q = (n_qatoms + BLOCK_SIZE - 1) / BLOCK_SIZE;
    int n_blocks_p = (n_patoms + BLOCK_SIZE - 1) / BLOCK_SIZE;
    if (!is_initialized) {
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

        is_initialized = true;
    }

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
    ctx.sync_all_to_device();

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