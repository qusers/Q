#include <vector>
#include <iostream>

#include "cuda/include/CudaContext.cuh"
#include "cuda/include/CudaNonbondedForce.cuh"
#include "cuda/include/CudaNonbondedQWForce.cuh"
#include "utils.h"

namespace CudaNonbondedQWForce {
bool is_initialized = false;

}  // namespace CudaNonbondedQWForce

void calc_nonbonded_qw_forces_host_v2() {
    using namespace CudaNonbondedQWForce;
    int nx = CudaContext::instance().h_q_atoms_list.size();
    int ny = CudaContext::instance().h_w_atoms_list.size();

    for (int state = 0; state < n_lambdas; state++) {
        auto result = calc_nonbonded_force_host(
            nx,
            ny,
            CudaContext::instance().d_q_atoms_list,
            CudaContext::instance().d_w_atoms_list,
            false,
            CudaContext::instance().d_q_charge_types + state * nx + nx * n_lambdas,
            CudaContext::instance().d_w_charge_types,
            CudaContext::instance().d_charge_table_all,
            CudaContext::instance().d_q_catype_types + state * nx + nx * n_lambdas,
            CudaContext::instance().d_w_catype_types,
            CudaContext::instance().d_catype_table_all);

        EQ_nonbond_qw[state].Uvdw += result.first / lambdas[state];
        EQ_nonbond_qw[state].Ucoul += result.second / lambdas[state];
        printf("Nonbonded QW Force State %d: Uvdw = %f, Ucoul = %f\n", state, result.first, result.second);
    }
}

void init_nonbonded_qw_force_kernel_data() {
    using namespace CudaNonbondedQWForce;
    if (is_initialized) return;
    is_initialized = true;
}

void cleanup_nonbonded_qw_force() {
    using namespace CudaNonbondedQWForce;
    if (!is_initialized) return;
    is_initialized = false;
}
