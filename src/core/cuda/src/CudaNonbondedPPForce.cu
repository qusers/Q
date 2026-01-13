#include <vector>

#include "cuda/include/CudaContext.cuh"
#include "cuda/include/CudaNonbondedForce.cuh"
#include "cuda/include/CudaNonbondedPPForce.cuh"
#include "utils.h"
#include <iostream>

namespace CudaNonbondedPPForce {
bool is_initialized = false;
int* d_x_idx_list = nullptr;
}  // namespace CudaNonbondedPPForce

void init_nonbonded_pp_force_kernel_data() {
    using namespace CudaNonbondedPPForce;

    if (!is_initialized) {
        is_initialized = true;
        std::vector<int> x_idx_list;
        x_idx_list.resize(n_patoms);
        for (int i = 0; i < n_patoms; i++) {
            x_idx_list[i] = p_atoms[i].a - 1;
        }

        check_cudaMalloc((void**)&d_x_idx_list, sizeof(int) * n_patoms);
        cudaMemcpy(
            d_x_idx_list,
            x_idx_list.data(),
            sizeof(int) * n_patoms,
            cudaMemcpyHostToDevice);
    }
}

void calc_nonbonded_pp_forces_host_v2() {
    int nx = n_patoms;
    int ny = n_patoms;

    auto result = calc_nonbonded_force_host(
        nx,
        ny,
        CudaNonbondedPPForce::d_x_idx_list,
        CudaNonbondedPPForce::d_x_idx_list,
        true  // symmetric
    );
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
    }
}