#include <vector>

#include "cuda/include/CudaContext.cuh"
#include "cuda/include/CudaNonbondedForce.cuh"
#include "cuda/include/CudaNonbondedQQForce.cuh"
#include "utils.h"

namespace CudaNonbondedQQForce {
bool is_initialized = false;
int* d_idx_list = nullptr;
int n = 0;

}  // namespace CudaNonbondedQQForce

void calc_nonbonded_qq_forces_host() {
    using namespace CudaNonbondedQQForce;

    for (int state = 0; state < n_lambdas; state++) {
        auto result = calc_nonbonded_force_host(
            n,
            n,
            d_idx_list,
            d_idx_list,
            true,
            CudaContext::instance().d_q_charge_types,
            CudaContext::instance().d_q_charge_types + state * n_qatoms,
            CudaContext::instance().d_charge_table_all,
            CudaContext::instance().d_q_catype_types,
            CudaContext::instance().d_q_catype_types + state * n_qatoms,
            CudaContext::instance().d_catype_table_all
        );

        EQ_nonbond_qq[state].Uvdw += result.first;
        EQ_nonbond_qq[state].Ucoul += result.second;
    }
}

void init_nonbonded_qq_force_kernel_data() {
    using namespace CudaNonbondedQQForce;
    if (is_initialized) return;

    std::vector<int> idx_list;
    for (int i = 0; i < n_qatoms; i++) {
        int id = q_atoms[i].a - 1;
        if (!excluded[id]) {
            idx_list.push_back(id);
        }
    }
    n = static_cast<int>(idx_list.size());

    check_cudaMalloc((void**)&d_idx_list, sizeof(int) * n);
    cudaMemcpy(d_idx_list, idx_list.data(), sizeof(int) * n, cudaMemcpyHostToDevice);

    is_initialized = true;
}

void cleanup_nonbonded_qq_force() {
    using namespace CudaNonbondedQQForce;
    if (!is_initialized) return;

    cudaFree(d_idx_list);

    d_idx_list = nullptr;
    n = 0;
    is_initialized = false;
}
