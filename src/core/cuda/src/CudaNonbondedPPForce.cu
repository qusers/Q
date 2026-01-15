#include <iostream>
#include <vector>

#include "cuda/include/CudaContext.cuh"
#include "cuda/include/CudaNonbondedForce.cuh"
#include "cuda/include/CudaNonbondedPPForce.cuh"
#include "utils.h"

namespace CudaNonbondedPPForce {
bool is_initialized = false;
int* d_x_idx_list = nullptr;
int nx;
}  // namespace CudaNonbondedPPForce

void init_nonbonded_pp_force_kernel_data() {
    using namespace CudaNonbondedPPForce;

    if (!is_initialized) {
        is_initialized = true;
        std::vector<int> x_idx_list;
        std::vector<int> x_charge_types;
        std::vector<int> x_catype_types;
        for (int i = 0; i < n_patoms; i++) {
            int id = p_atoms[i].a - 1;
            if (!excluded[id]) {
                x_idx_list.push_back(id);
            }
        }
        nx = x_idx_list.size();

        check_cudaMalloc((void**)&d_x_idx_list, sizeof(int) * x_idx_list.size());
        cudaMemcpy(
            d_x_idx_list,
            x_idx_list.data(),
            sizeof(int) * x_idx_list.size(),
            cudaMemcpyHostToDevice);
    }
}

void calc_nonbonded_pp_forces_host_v2() {
    int nx = CudaNonbondedPPForce::nx;
    int ny = nx;

    auto result = calc_nonbonded_force_host(
        nx,
        ny,
        CudaNonbondedPPForce::d_x_idx_list,
        CudaNonbondedPPForce::d_x_idx_list,
        true,  // symmetric
        CudaContext::instance().d_p_charge_types,
        CudaContext::instance().d_p_charge_types,
        CudaContext::instance().d_charge_table_all,
        CudaContext::instance().d_p_catype_types,
        CudaContext::instance().d_p_catype_types,
        CudaContext::instance().d_catype_table_all);
    printf("Nonbonded PP Force: Uvdw = %f, Ucoul = %f\n", result.first, result.second);
    E_nonbond_pp.Uvdw += result.first;
    E_nonbond_pp.Ucoul += result.second;
}

void cleanup_nonbonded_pp_force() {
    using namespace CudaNonbondedPPForce;

    if (is_initialized) {
        is_initialized = false;
        cudaFree(d_x_idx_list);
        d_x_idx_list = nullptr;
        nx = 0;
    }
}
