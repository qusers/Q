#include <vector>

#include "cuda/include/CudaContext.cuh"
#include "cuda/include/CudaNonbondedForce.cuh"
#include "cuda/include/CudaNonbondedPWForce.cuh"
#include <iostream>

namespace CudaNonbondedPWForce {
// Declare any necessary static variables or device pointers here
bool is_initialized = false;
int* d_x_idx_list = nullptr;
int* d_y_idx_list = nullptr;

}  // namespace CudaNonbondedPWForce
void calc_nonbonded_pw_forces_host_v2() {
    int nx = n_patoms;
    int ny = n_waters * 3;

    using namespace CudaNonbondedPWForce;
    auto result = calc_nonbonded_force_host(nx, ny, d_x_idx_list, d_y_idx_list, false);
    printf("Nonbonded PW Force (Host) - VdW: %f, Coulomb: %f\n", result.first, result.second);

    E_nonbond_pw.Uvdw = result.first;
    E_nonbond_pw.Ucoul = result.second;
}

void init_nonbonded_pw_force_kernel_data() {
    using namespace CudaNonbondedPWForce;
    if (!is_initialized) {
        std::vector<int> x_idx_list(n_patoms);
        std::vector<int> y_idx_list(n_waters * 3);

        for (int i = 0; i < n_patoms; i++) {
            x_idx_list[i] = p_atoms[i].a - 1;
        }

        for (int i = n_atoms_solute; i < n_atoms; i++) {
            y_idx_list[i - n_atoms_solute] = i;
        }
        check_cudaMalloc((void**)&d_x_idx_list, sizeof(int) * n_patoms);
        check_cudaMalloc((void**)&d_y_idx_list, sizeof(int) * n_waters * 3);
        cudaMemcpy(d_x_idx_list, x_idx_list.data(), sizeof(int) * n_patoms, cudaMemcpyHostToDevice);
        cudaMemcpy(d_y_idx_list, y_idx_list.data(), sizeof(int) * n_waters * 3, cudaMemcpyHostToDevice);
        is_initialized = true;
    }
}

void cleanup_nonbonded_pw_force() {
    using namespace CudaNonbondedPWForce;
    if (is_initialized) {
        cudaFree(d_x_idx_list);
        cudaFree(d_y_idx_list);
        d_x_idx_list = nullptr;
        d_y_idx_list = nullptr;
        is_initialized = false;
    }
}