#include <iostream>
#include <vector>

#include "common/include/context.h"
#include "cuda/include/cuda_context.cuh"
#include "cuda/include/cuda_nonbonded_force.cuh"
#include "cuda/include/cuda_nonbonded_ww_force.cuh"
namespace CudaNonbondedWWForce {

bool is_initialized = false;
}  // namespace CudaNonbondedWWForce
void calc_nonbonded_ww_forces_host_v2() {
    using namespace CudaNonbondedWWForce;

    int nx = CudaContext::instance().h_w_atoms_list.size();
    int ny = nx;
    auto result = calc_nonbonded_force_host(
        nx, ny,
        CudaContext::instance().d_w_atoms_list,
        CudaContext::instance().d_w_atoms_list,
        true,
        CudaContext::instance().d_w_charge_types,
        CudaContext::instance().d_w_charge_types,
        CudaContext::instance().d_w_catype_types,
        CudaContext::instance().d_w_catype_types,
        true);
    Context& host = Context::instance();
    host.E_nonbond_ww.Uvdw = result.first;
    host.E_nonbond_ww.Ucoul = result.second;
}

void init_nonbonded_ww_force_kernel_data() {
    using namespace CudaNonbondedWWForce;
    if (!is_initialized) {
        is_initialized = true;
    }
}

void cleanup_nonbonded_ww_force() {
    using namespace CudaNonbondedWWForce;
    if (is_initialized) {
        is_initialized = false;
    }
}
