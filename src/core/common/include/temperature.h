#pragma once

#include "context.h"
#include "shake.h"

enum TResult { R_TEMP,
               R_TFREE,
               R_TSCALE_SOL,
               R_TSCALE_SLV,
               R_UKIN,
               N_TRESULT };
               
struct TemperatureData {
    // init
    int Ndegf, Ndegfree, Ndegf_solute, Ndegf_solvent, Ndegfree_solute, Ndegfree_solvent;
    double tau_T = 0;
    bool separate_scaling = false;

    // calculate, every step it will change.
    // device slot read by leapfrog kernel, host slot for logging
    std::unique_ptr<HostDeviceBuffer<double>> results;
};

class Temperature {
   public:
    virtual ~Temperature() = default;
    void init(Context& ctx, const Shake& shake);
    virtual void calc(Context& ctx) = 0;
    virtual void sync_for_output(Context& ctx, int replica = 0);
    virtual void cleanup() {}
    const TemperatureData& data() const { return data_; }

   protected:
    TemperatureData data_;
    virtual void init_backend(Context& ctx) {}
};