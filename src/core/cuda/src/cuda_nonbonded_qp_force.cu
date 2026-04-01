#include <iostream>
#include <map>
#include <vector>

#include "common/include/context.h"
#include "cuda/include/cuda_context.cuh"
#include "cuda/include/cuda_nonbonded_force.cuh"
#include "cuda/include/cuda_nonbonded_qp_force.cuh"
#include "cuda_utility.cuh"

namespace CudaNonbondedQPForce {
bool is_initialized = false;

}  // namespace CudaNonbondedQPForce

void calc_nonbonded_qp_forces_host_v2() {
    using namespace CudaNonbondedQPForce;
    int nx = CudaContext::instance().h_q_atoms_list.size();
    int ny = CudaContext::instance().h_p_atoms_list.size();
    
    Context& host = Context::instance();
    for (int state = 0; state < host.n_lambdas; state++) {
        auto result = calc_nonbonded_force_host(
            nx,
            ny,
            CudaContext::instance().d_q_atoms_list,
            CudaContext::instance().d_p_atoms_list,
            false,
            CudaContext::instance().d_q_charge_types + state * nx,
            CudaContext::instance().d_p_charge_types,
            CudaContext::instance().d_charge_table_all,
            CudaContext::instance().d_q_catype_types + state * nx,
            CudaContext::instance().d_p_catype_types,
            CudaContext::instance().d_catype_table_all, false, host.lambdas[state]);

        host.EQ_nonbond_qp[state].Uvdw = result.first / host.lambdas[state];
        host.EQ_nonbond_qp[state].Ucoul = result.second / host.lambdas[state];
        // printf("Nonbonded QP Force State %d: Uvdw = %f, Ucoul = %f\n", state, result.first, result.second);
    }
}

void init_nonbonded_qp_force_kernel_data() {
    using namespace CudaNonbondedQPForce;
    if (!is_initialized) {
        is_initialized = true;
    }
}

void cleanup_nonbonded_qp_force() {
    using namespace CudaNonbondedQPForce;
    if (is_initialized) {
        is_initialized = false;
    }
}
