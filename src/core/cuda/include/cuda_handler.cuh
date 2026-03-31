#pragma once

#include "handler.h"
#include "cuda_context.cuh"

class CudaHandler : public Handler {
   public:
    static CudaHandler& instance() {
        static CudaHandler handler;
        return handler;
    }

    // Allocate device resources and prepare kernel state.
    void initialize() override;

    // Release device resources.
    void shutdown() override;


   protected:
    bool initialized_ = false;
    CudaContext& ctx_ = CudaContext::instance();

    void calc_internal_forces(int iteration) override;
    void calc_nonbonded_forces() override;
    void calc_temperature() override;
    void calc_leapfrog() override;

    void reset_energies() override;
};
