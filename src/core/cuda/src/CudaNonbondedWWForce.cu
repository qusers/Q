#include <iostream>
#include <vector>

#include "cuda/include/CudaContext.cuh"
#include "cuda/include/CudaNonbondedForce.cuh"
#include "cuda/include/CudaNonbondedWWForce.cuh"
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
        CudaContext::instance().d_charge_table_all,
        CudaContext::instance().d_w_catype_types,
        CudaContext::instance().d_w_catype_types,
        CudaContext::instance().d_catype_table_all);
    E_nonbond_ww.Uvdw = result.first;
    E_nonbond_ww.Ucoul = result.second;
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
