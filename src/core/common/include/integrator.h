#pragma once

#include "constraint_force.h"
#include "context.h"
#include "temperature.h"

struct IntegratorData {
    std::unique_ptr<HostDeviceBuffer<coord_t>> xcoords;
};

class Integrator {
   public:
    virtual ~Integrator() = default;

    // Capture collaborators + backend setup. constraint force is applied inside every
    // step; temperature supplies the per-step Tscale results.
    void init(Context& ctx, ConstraintForce& constraint_force, const Temperature& temperature);

    virtual void step(Context& ctx) = 0;  // one leapfrog step
    virtual void cleanup() {}

    const IntegratorData& data() const { return data_; }

   protected:
    IntegratorData data_;
    ConstraintForce* constraint_force_ = nullptr;  // not owned
    const Temperature* temperature_ = nullptr;     // not owned; source of Tscale

    virtual void init_backend(Context& ctx) {}
};