#pragma once

#include "context.h"
#include "shake.h"
#include "temperature.h"

struct IntegratorData {
};

class Integrator {
   public:
    virtual ~Integrator() = default;

    // Capture collaborators + backend setup. shake is applied inside every
    // step; temperature supplies the per-step Tscale results.
    void init(Context& ctx, Shake& shake, const Temperature& temperature);

    virtual void step(Context& ctx) = 0;  // one leapfrog step
    virtual void cleanup() {}

    const IntegratorData& data() const { return data_; }

   protected:
    IntegratorData data_;
    Shake* shake_ = nullptr;                    // not owned
    const Temperature* temperature_ = nullptr;  // not owned; source of Tscale

    virtual void init_backend(Context& ctx) {}
};