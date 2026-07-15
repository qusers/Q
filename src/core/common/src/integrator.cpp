#include "integrator.h"

void Integrator::init(Context& ctx, Shake& shake, const Temperature& temperature) {
    shake_ = &shake;
    temperature_ = &temperature;
    data_.xcoords = std::make_unique<HostDeviceBuffer<coord_t>>(ctx.n_atoms, true, ctx.command_info.requested_gpu);
    init_backend(ctx);
}
