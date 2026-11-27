#include "integrator.h"

void Integrator::init(Context& ctx, ConstraintForce& constraint_force, const Temperature& temperature) {
    constraint_force_ = &constraint_force;
    temperature_ = &temperature;
    data_.xcoords = std::make_unique<HostDeviceBuffer<coord_t>>(ctx.n_atoms, true, ctx.command_info.requested_gpu);
    init_backend(ctx);
}
