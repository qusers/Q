#include <iostream>
#include <vector>

#include "common/include/context.h"
#include "cuda/include/cuda_nonbonded_force.cuh"
#include "cuda/include/cuda_nonbonded_pw_force.cuh"

namespace CudaNonbondedPWForce {
// Declare any necessary static variables or device pointers here
bool is_initialized = false;

}  // namespace CudaNonbondedPWForce
void calc_nonbonded_pw_forces_host_v2() {
    using namespace CudaNonbondedPWForce;

    Context& host = Context::instance();
    int nx = static_cast<int>(host.p_atoms_list->length);
    int ny = static_cast<int>(host.w_atoms_list->length);
    auto result = calc_nonbonded_force_host(
        nx, ny,
        host.p_atoms_list->gpu_data_p,
        host.w_atoms_list->gpu_data_p,
        false,
        host.p_charge_types->gpu_data_p,
        host.w_charge_types->gpu_data_p,
        host.p_catype_types->gpu_data_p,
        host.w_catype_types->gpu_data_p,
        false);
    // printf("Nonbonded PW Force (Host) - VdW: %f, Coulomb: %f\n", result.first, result.second);
    host.E_nonbond_pw.Uvdw = result.first;
    host.E_nonbond_pw.Ucoul = result.second;
}

void init_nonbonded_pw_force_kernel_data() {
    using namespace CudaNonbondedPWForce;
    if (!is_initialized) {
        is_initialized = true;
    }
}

void cleanup_nonbonded_pw_force() {
    using namespace CudaNonbondedPWForce;
    if (is_initialized) {
        is_initialized = false;
    }
}
