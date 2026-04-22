#include <vector>
#include <iostream>

#include "common/include/context.h"
#include "cuda/include/cuda_nonbonded_force.cuh"
#include "cuda/include/cuda_nonbonded_qw_force.cuh"
#include "cuda_utility.cuh"

namespace CudaNonbondedQWForce {
bool is_initialized = false;

}  // namespace CudaNonbondedQWForce

void calc_nonbonded_qw_forces_host_v2() {
    using namespace CudaNonbondedQWForce;

    Context& host = Context::instance();
    int nx = static_cast<int>(host.q_atoms_list->length);
    int ny = static_cast<int>(host.w_atoms_list->length);
    auto *lambdas = host.lambdas->cpu_data_p;
    for (int state = 0; state < host.n_lambdas; state++) {
        auto result = calc_nonbonded_force_host(
            nx,
            ny,
            host.q_atoms_list->gpu_data_p,
            host.w_atoms_list->gpu_data_p,
            false,
            host.q_charge_types->gpu_data_p + state * nx,
            host.w_charge_types->gpu_data_p,
            host.q_catype_types->gpu_data_p + state * nx,
            host.w_catype_types->gpu_data_p,
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
