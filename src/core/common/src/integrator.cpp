#include "integrator.h"

void Integrator::init(Context& ctx, Shake& shake, const Temperature& temperature) {
    shake_ = &shake;
    temperature_ = &temperature;
    init_backend(ctx);
}

