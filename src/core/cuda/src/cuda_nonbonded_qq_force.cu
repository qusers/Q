#include <vector>
#include <iostream>

#include "common/include/context.h"
#include "cuda/include/cuda_context.cuh"
#include "cuda/include/cuda_nonbonded_force.cuh"
#include "cuda/include/cuda_nonbonded_qq_force.cuh"
#include "cuda_utility.cuh"

namespace CudaNonbondedQQForce {
bool is_initialized = false;

}  // namespace CudaNonbondedQQForce

void calc_nonbonded_qq_forces_host() {
    using namespace CudaNonbondedQQForce;

    Context& host = Context::instance();
    int n = CudaContext::instance().h_q_atoms_list.size();
    for (int state = 0; state < host.n_lambdas; state++) {
        auto result = calc_nonbonded_force_host(
            n,
            n,
            CudaContext::instance().d_q_atoms_list,
            CudaContext::instance().d_q_atoms_list,
            true,
            CudaContext::instance().d_q_charge_types + state * n,
            CudaContext::instance().d_q_charge_types + state * n,
            CudaContext::instance().d_q_catype_types + state * n,
            CudaContext::instance().d_q_catype_types + state * n,
            false, host.lambdas[state]);

        host.EQ_nonbond_qq[state].Uvdw = result.first / host.lambdas[state];
        host.EQ_nonbond_qq[state].Ucoul = result.second / host.lambdas[state];
    }
}

void init_nonbonded_qq_force_kernel_data() {
    using namespace CudaNonbondedQQForce;
    if (is_initialized) return;
    is_initialized = true;
}

void cleanup_nonbonded_qq_force() {
    using namespace CudaNonbondedQQForce;
    if (!is_initialized) return;
    is_initialized = false;
}
