#include "cuda_shake_v2.cuh"

void CudaShakeV2::apply(Context& ctx) {
    ctx.coords->download();
    ctx.xcoords->download();
    CpuShake::apply(ctx);
    ctx.coords->upload();
}