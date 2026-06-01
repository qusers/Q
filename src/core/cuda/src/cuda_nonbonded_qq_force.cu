#include <vector>
#include <iostream>

#include "common/include/context.h"
#include "cuda/include/cuda_nonbonded_force.cuh"
#include "cuda/include/cuda_nonbonded_qq_force.cuh"
#include "cuda_utility.cuh"

namespace CudaNonbondedQQForce {
bool is_initialized = false;

}  // namespace CudaNonbondedQQForce

void calc_nonbonded_qq_forces_host() {
    using namespace CudaNonbondedQQForce;

    Context& host = Context::instance();
    auto *lambdas = host.lambdas->cpu_data_p;
    int n = static_cast<int>(host.q_atoms_list->length);
    for (int state = 0; state < host.n_lambdas; state++) {
        auto result = calc_nonbonded_force_host(
            n,
            n,
            host.q_atoms_list->gpu_data_p,
            host.q_atoms_list->gpu_data_p,
            true,
            host.q_charge_types->gpu_data_p + state * n,
            host.q_charge_types->gpu_data_p + state * n,
            host.q_catype_types->gpu_data_p + state * n,
            host.q_catype_types->gpu_data_p + state * n,
            false, lambdas[state]);

        host.EQ_nonbond_qq[state].Uvdw = result.first / lambdas[state];
        host.EQ_nonbond_qq[state].Ucoul = result.second / lambdas[state];
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
