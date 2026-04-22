#include <iostream>
#include <vector>

#include "common/include/context.h"
#include "cuda/include/cuda_nonbonded_force.cuh"
#include "cuda/include/cuda_nonbonded_ww_force.cuh"
namespace CudaNonbondedWWForce {

bool is_initialized = false;
}  // namespace CudaNonbondedWWForce
void calc_nonbonded_ww_forces_host_v2() {
    using namespace CudaNonbondedWWForce;

    Context& host = Context::instance();
    int nx = static_cast<int>(host.w_atoms_list->length);
    int ny = nx;
    auto result = calc_nonbonded_force_host(
        nx, ny,
        host.w_atoms_list->gpu_data_p,
        host.w_atoms_list->gpu_data_p,
        true,
        host.w_charge_types->gpu_data_p,
        host.w_charge_types->gpu_data_p,
        host.w_catype_types->gpu_data_p,
        host.w_catype_types->gpu_data_p,
        true);
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
