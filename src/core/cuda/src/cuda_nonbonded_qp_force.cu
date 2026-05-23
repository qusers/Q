#include <iostream>
#include <map>
#include <vector>

#include "common/include/context.h"
#include "cuda/include/cuda_nonbonded_force.cuh"
#include "cuda/include/cuda_nonbonded_qp_force.cuh"
#include "cuda_utility.cuh"

namespace CudaNonbondedQPForce {
bool is_initialized = false;

}  // namespace CudaNonbondedQPForce

void calc_nonbonded_qp_forces_host_v2() {
    using namespace CudaNonbondedQPForce;
    Context& host = Context::instance();
    int nx = static_cast<int>(host.q_atoms_list->length);
    int ny = static_cast<int>(host.p_atoms_list->length);
    auto *lambdas = host.lambdas->cpu_data_p;
    for (int state = 0; state < host.n_lambdas; state++) {
        const bool softcore_enabled = host.n_qsoftcores > 0;
        const real_t* softcore_values = softcore_enabled
            ? host.q_softcore_values->gpu_data_p + state * host.n_qatoms
            : nullptr;
        auto result = calc_nonbonded_force_host(
            nx,
            ny,
            host.q_atoms_list->gpu_data_p,
            host.p_atoms_list->gpu_data_p,
            false,
            host.q_charge_types->gpu_data_p + state * nx,
            host.p_charge_types->gpu_data_p,
            host.q_catype_types->gpu_data_p + state * nx,
            host.p_catype_types->gpu_data_p,
            false,
            lambdas[state],
            softcore_enabled,
            softcore_values,
            host.atom_to_qi_lookup->gpu_data_p,
            host.q_softcore_use_max_potential);

        host.EQ_nonbond_qp[state].Uvdw = result.first / lambdas[state];
        host.EQ_nonbond_qp[state].Ucoul = result.second / lambdas[state];
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
