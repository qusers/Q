#include <vector>
#include <iostream>

#include "cuda/include/CudaContext.cuh"
#include "cuda/include/CudaNonbondedForce.cuh"
#include "cuda/include/CudaNonbondedQQForce.cuh"
#include "utils.h"

namespace CudaNonbondedQQForce {
bool is_initialized = false;

}  // namespace CudaNonbondedQQForce

void calc_nonbonded_qq_forces_host() {
    using namespace CudaNonbondedQQForce;

    int n = CudaContext::instance().h_q_atoms_list.size();
    for (int state = 0; state < n_lambdas; state++) {
        auto result = calc_nonbonded_force_host(
            n,
            n,
            CudaContext::instance().d_q_atoms_list,
            CudaContext::instance().d_q_atoms_list,
            true,
            CudaContext::instance().d_q_charge_types + state * n,
            CudaContext::instance().d_q_charge_types + state * n + n * n_lambdas,
            CudaContext::instance().d_charge_table_all,
            CudaContext::instance().d_q_catype_types + state * n,
            CudaContext::instance().d_q_catype_types + state * n + n * n_lambdas,
            CudaContext::instance().d_catype_table_all);

        EQ_nonbond_qq[state].Uvdw += result.first / lambdas[state];
        EQ_nonbond_qq[state].Ucoul += result.second / lambdas[state];
        printf("Nonbonded QQ Force State %d: Uvdw = %f, Ucoul = %f\n", state, result.first, result.second);
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
