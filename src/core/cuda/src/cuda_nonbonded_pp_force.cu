#include <iostream>
#include <vector>

#include "common/include/context.h"
#include "cuda/include/cuda_nonbonded_force.cuh"
#include "cuda/include/cuda_nonbonded_pp_force.cuh"
#include "cuda_utility.cuh"

namespace CudaNonbondedPPForce {
bool is_initialized = false;
}  // namespace CudaNonbondedPPForce

void init_nonbonded_pp_force_kernel_data() {
    using namespace CudaNonbondedPPForce;

    if (!is_initialized) {
        is_initialized = true;
    }
}

void calc_nonbonded_pp_forces_host_v2() {
    Context& host = Context::instance();
    int nx = static_cast<int>(host.p_atoms_list->length);
    int ny = nx;

    auto result = calc_nonbonded_force_host(
        nx,
        ny,
        host.p_atoms_list->gpu_data_p,
        host.p_atoms_list->gpu_data_p,
        true,  // symmetric
        host.p_charge_types->gpu_data_p,
        host.p_charge_types->gpu_data_p,
        host.p_catype_types->gpu_data_p,
        host.p_catype_types->gpu_data_p,
        false);
    host.E_nonbond_pp.Uvdw = result.first;
    host.E_nonbond_pp.Ucoul = result.second;
}

void cleanup_nonbonded_pp_force() {
    using namespace CudaNonbondedPPForce;

    if (is_initialized) {
        is_initialized = false;
    }
}
