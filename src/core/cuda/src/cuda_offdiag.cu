#include "offdiag.h"

namespace {

__global__ void offdiag_kernel(
    int n,
    const offdiag_term_t* terms,
    const coord_t* coords,
    offdiag_record_t* records) {
    const int i = blockIdx.x * blockDim.x + threadIdx.x;
    if (i >= n) return;

    records[i] = calc_offdiag(terms[i], coords);
}

}  // namespace

void calc_offdiags_cuda(Context& ctx) {
    const int n = ctx.n_offdiags();
    if (n == 0) return;

    const int block = 128;
    const int grid = (n + block - 1) / block;

    offdiag_kernel<<<grid, block>>>(
        n,
        ctx.offdiag_terms->gpu_data_p,
        ctx.coords->gpu_data_p,
        ctx.offdiag_records->gpu_data_p);

    check_cuda(cudaGetLastError());
}