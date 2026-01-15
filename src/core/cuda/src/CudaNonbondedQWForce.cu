#include <vector>

#include "cuda/include/CudaContext.cuh"
#include "cuda/include/CudaNonbondedForce.cuh"
#include "cuda/include/CudaNonbondedQWForce.cuh"
#include "utils.h"

namespace CudaNonbondedQWForce {
bool is_initialized = false;
int* d_x_idx_list = nullptr;  // q atoms
int* d_y_idx_list = nullptr;  // water atoms
int nx = 0, ny = 0;

}  // namespace CudaNonbondedQWForce

void calc_nonbonded_qw_forces_host_v2() {
    using namespace CudaNonbondedQWForce;

    for (int state = 0; state < n_lambdas; state++) {
        auto result = calc_nonbonded_force_host(
            nx,
            ny,
            d_x_idx_list,
            d_y_idx_list,
            false,
            CudaContext::instance().d_q_charge_types + state * n_qatoms,
            CudaContext::instance().d_w_charge_types,
            CudaContext::instance().d_charge_table_all,
            CudaContext::instance().d_q_catype_types + state * n_qatoms,
            CudaContext::instance().d_w_catype_types,
            CudaContext::instance().d_catype_table_all);

        EQ_nonbond_qw[state].Uvdw += result.first;
        EQ_nonbond_qw[state].Ucoul += result.second;
    }
}

void init_nonbonded_qw_force_kernel_data() {
    using namespace CudaNonbondedQWForce;
    if (is_initialized) return;

    std::vector<int> x_idx_list;
    std::vector<int> y_idx_list;

    // q atoms
    for (int i = 0; i < n_qatoms; i++) {
        int id = q_atoms[i].a - 1;
        if (!excluded[id]) {
            x_idx_list.push_back(id);
        }
    }
    // water atoms
    for (int i = n_atoms_solute; i < n_atoms; i++) {
        if (!excluded[i]) {
            y_idx_list.push_back(i);
        }
    }

    nx = static_cast<int>(x_idx_list.size());
    ny = static_cast<int>(y_idx_list.size());

    check_cudaMalloc((void**)&d_x_idx_list, sizeof(int) * nx);
    check_cudaMalloc((void**)&d_y_idx_list, sizeof(int) * ny);
    cudaMemcpy(d_x_idx_list, x_idx_list.data(), sizeof(int) * nx, cudaMemcpyHostToDevice);
    cudaMemcpy(d_y_idx_list, y_idx_list.data(), sizeof(int) * ny, cudaMemcpyHostToDevice);

    is_initialized = true;
}

void cleanup_nonbonded_qw_force() {
    using namespace CudaNonbondedQWForce;
    if (!is_initialized) return;

    cudaFree(d_x_idx_list);
    cudaFree(d_y_idx_list);

    d_x_idx_list = nullptr;
    d_y_idx_list = nullptr;
    nx = ny = 0;
    is_initialized = false;
}
