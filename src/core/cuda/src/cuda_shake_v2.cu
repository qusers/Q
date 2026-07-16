#include "cuda_shake_v2.cuh"
void CudaShakeV2::apply(Context& ctx, HostDeviceBuffer<coord_t>& xcoords) {
    ctx.coords->download();
    CpuShake::apply(ctx, xcoords);
    ctx.coords->upload();
}