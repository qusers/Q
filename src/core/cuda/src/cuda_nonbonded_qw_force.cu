#include <vector>
#include <iostream>

#include "common/include/context.h"
#include "cuda/include/cuda_context.cuh"
#include "cuda/include/cuda_nonbonded_force.cuh"
#include "cuda/include/cuda_nonbonded_qw_force.cuh"
#include "cuda_utility.cuh"

namespace CudaNonbondedQWForce {
bool is_initialized = false;

}  // namespace CudaNonbondedQWForce

void calc_nonbonded_qw_forces_host_v2() {
    using namespace CudaNonbondedQWForce;
    int nx = CudaContext::instance().h_q_atoms_list.size();
    int ny = CudaContext::instance().h_w_atoms_list.size();

    Context& host = Context::instance();
    auto *lambdas = host.lambdas->cpu_data_p;
    for (int state = 0; state < host.n_lambdas; state++) {
        auto result = calc_nonbonded_force_host(
            nx,
            ny,
            CudaContext::instance().d_q_atoms_list,
            CudaContext::instance().d_w_atoms_list,
            false,
            CudaContext::instance().d_q_charge_types + state * nx,
            CudaContext::instance().d_w_charge_types,
            CudaContext::instance().d_q_catype_types + state * nx,
            CudaContext::instance().d_w_catype_types,
            true, lambdas[state]);

        host.EQ_nonbond_qw[state].Uvdw = result.first / lambdas[state];
        host.EQ_nonbond_qw[state].Ucoul = result.second / lambdas[state];
        // printf("Nonbonded QW Force State %d: Uvdw = %f, Ucoul = %f\n", state, result.first, result.second);
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
