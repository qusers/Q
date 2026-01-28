#pragma once

#include "CudaContext.cuh"

class CudaHandler {
   public:
    static CudaHandler& instance() {
        static CudaHandler handler;
        return handler;
    }

    // Allocate device resources and prepare kernel state.
    void initialize();

    // Release device resources.
    void shutdown();

    // Run a single GPU simulation iteration (mirrors calc_integration_step).
    void run_iteration(int iteration);

    void run(int num_iterations);

   private:
    CudaHandler() = default;
    CudaHandler(const CudaHandler&) = delete;
    CudaHandler& operator=(const CudaHandler&) = delete;

    bool initialized_ = false;
    CudaContext& ctx_ = CudaContext::instance();

    void calc_internal_forces(int iteration);
    void calc_nonbonded_forces();
    void update_energy_totals();
    void print_outputs(int iteration);
    void reset_energies();
};
