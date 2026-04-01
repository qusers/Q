#include <iostream>
#include <vector>

#include "common/include/context.h"
#include "cuda/include/cuda_context.cuh"
#include "cuda/include/cuda_nonbonded_force.cuh"
#include "cuda/include/cuda_nonbonded_pw_force.cuh"

namespace CudaNonbondedPWForce {
// Declare any necessary static variables or device pointers here
bool is_initialized = false;

}  // namespace CudaNonbondedPWForce
void calc_nonbonded_pw_forces_host_v2() {
    using namespace CudaNonbondedPWForce;

    int nx = CudaContext::instance().h_p_atoms_list.size();
    int ny = CudaContext::instance().h_w_atoms_list.size();
    auto result = calc_nonbonded_force_host(
        nx, ny,
        CudaContext::instance().d_p_atoms_list,
        CudaContext::instance().d_w_atoms_list,
        false,
        CudaContext::instance().d_p_charge_types,
        CudaContext::instance().d_w_charge_types,
        CudaContext::instance().d_charge_table_all,
        CudaContext::instance().d_p_catype_types,
        CudaContext::instance().d_w_catype_types,
        CudaContext::instance().d_catype_table_all, false);
    // printf("Nonbonded PW Force (Host) - VdW: %f, Coulomb: %f\n", result.first, result.second);

    Context& host = Context::instance();
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
