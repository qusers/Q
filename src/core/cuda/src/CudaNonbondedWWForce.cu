#include <iostream>
#include <vector>

#include "cuda/include/CudaContext.cuh"
#include "cuda/include/CudaNonbondedForce.cuh"
#include "cuda/include/CudaNonbondedWWForce.cuh"
namespace CudaNonbondedWWForce {

bool is_initialized = false;
int* d_x_idx_list = nullptr;
}  // namespace CudaNonbondedWWForce
void calc_nonbonded_ww_forces_host_v2() {
    using namespace CudaNonbondedWWForce;
    auto result = calc_nonbonded_force_host(n_waters * 3, n_waters * 3, d_x_idx_list, d_x_idx_list, true);
    E_nonbond_ww.Uvdw = result.first;
    E_nonbond_ww.Ucoul = result.second;
}

void init_nonbonded_ww_force_kernel_data() {
    using namespace CudaNonbondedWWForce;
    if (!is_initialized) {
        if (crg_ow == 0 || crg_hw == 0) {
            // Initialize water charges once (used by QW kernels too)
            // todo: Don't do that in here ...
            ccharge_t ccharge_ow = ccharges[charges[n_atoms_solute].code - 1];
            ccharge_t ccharge_hw = ccharges[charges[n_atoms_solute + 1].code - 1];
            crg_ow = ccharge_ow.charge;
            crg_hw = ccharge_hw.charge;
        }

        std::vector<int> x_idx_list(n_waters * 3);
        for (int i = n_atoms_solute; i < n_atoms; i++) {
            x_idx_list[i - n_atoms_solute] = i;
        }
        check_cudaMalloc((void**)&d_x_idx_list, sizeof(int) * n_waters * 3);
        cudaMemcpy(d_x_idx_list, x_idx_list.data(), sizeof(int) * n_waters * 3, cudaMemcpyHostToDevice);
        is_initialized = true;
    }
}

void cleanup_nonbonded_ww_force() {
    using namespace CudaNonbondedWWForce;
    if (is_initialized) {
        cudaFree(d_x_idx_list);
        d_x_idx_list = nullptr;
        is_initialized = false;
    }
}
