#include <iostream>
#include <map>
#include <vector>

#include "cuda/include/CudaContext.cuh"
#include "cuda/include/CudaNonbondedQPForce.cuh"
#include "cuda/include/CudaNonbondedForce.cuh"
#include "system.h"
#include "utils.h"

namespace CudaNonbondedQPForce {
bool is_initialized = false;
int* d_x_idx_list = nullptr;
int* d_y_idx_list = nullptr;

int nx, ny;

}  // namespace CudaNonbondedQPForce

void calc_nonbonded_qp_forces_host_v2() {
    using namespace CudaNonbondedQPForce;
    for (int state = 0; state < n_lambdas; state++) {
        auto result = calc_nonbonded_force_host(
            nx,
            ny,
            d_x_idx_list,
            d_y_idx_list,
            false,
            CudaContext::instance().d_q_charge_types + state * n_qatoms,
            CudaContext::instance().d_p_charge_types,
            CudaContext::instance().d_charge_table_all,
            CudaContext::instance().d_q_catype_types + state * n_qatoms,
            CudaContext::instance().d_p_catype_types,
            CudaContext::instance().d_catype_table_all);

        EQ_nonbond_qp[state].Uvdw += result.first;
        EQ_nonbond_qp[state].Ucoul += result.second;
    }
}

void init_nonbonded_qp_force_kernel_data() {
    using namespace CudaNonbondedQPForce;
    if (!is_initialized) {
        is_initialized = true;

        std::vector<int> x_idx_list;
        std::vector<int> y_idx_list;

        for (int i = 0; i < n_qatoms; i++) {
            int id = q_atoms[i].a - 1;
            if (!excluded[id]) {
                x_idx_list.push_back(id);
            }
        }

        for (int i = 0; i < n_patoms; i++) {
            int id = p_atoms[i].a - 1;
            if (!excluded[id]) {
                y_idx_list.push_back(id);
            }
        }

        nx = static_cast<int>(x_idx_list.size());
        ny = static_cast<int>(y_idx_list.size());

        check_cudaMalloc((void**)&d_x_idx_list, sizeof(int) * nx);
        check_cudaMalloc((void**)&d_y_idx_list, sizeof(int) * ny);
        cudaMemcpy(d_x_idx_list, x_idx_list.data(), sizeof(int) * nx, cudaMemcpyHostToDevice);
        cudaMemcpy(d_y_idx_list, y_idx_list.data(), sizeof(int) * ny, cudaMemcpyHostToDevice);
    }
}

void cleanup_nonbonded_qp_force() {
    using namespace CudaNonbondedQPForce;
    if (is_initialized) {
        cudaFree(d_x_idx_list);
        cudaFree(d_y_idx_list);
        d_x_idx_list = nullptr;
        d_y_idx_list = nullptr;
        nx = 0;
        ny = 0;
        is_initialized = false;
    }
}
