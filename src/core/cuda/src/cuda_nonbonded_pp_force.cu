#include <iostream>
#include <vector>

#include "common/include/context.h"
#include "cuda/include/cuda_context.cuh"
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
    int nx = CudaContext::instance().h_p_atoms_list.size();
    int ny = nx;

    auto result = calc_nonbonded_force_host(
        nx,
        ny,
        CudaContext::instance().d_p_atoms_list,
        CudaContext::instance().d_p_atoms_list,
        true,  // symmetric
        CudaContext::instance().d_p_charge_types,
        CudaContext::instance().d_p_charge_types,
        CudaContext::instance().d_charge_table_all,
        CudaContext::instance().d_p_catype_types,
        CudaContext::instance().d_p_catype_types,
        CudaContext::instance().d_catype_table_all, false);

    Context& host = Context::instance();
    host.E_nonbond_pp.Uvdw = result.first;
    host.E_nonbond_pp.Ucoul = result.second;
}

void cleanup_nonbonded_pp_force() {
    using namespace CudaNonbondedPPForce;

    if (is_initialized) {
        is_initialized = false;
    }
}
