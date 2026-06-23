#pragma once

#include "handler.h"

class CudaHandler : public Handler {
   public:
    static CudaHandler& instance() {
        static CudaHandler handler;
        return handler;
    }

    // Release device resources.
    void shutdown() override;


   protected:
    bool initialized_ = false;

    void initialize_backend() override;
    void calc_internal_forces(int iteration) override;
    void calc_nonbonded_forces() override;
    void calc_temperature() override;
    void calc_leapfrog() override;

    void reset_energies() override;
    std::unique_ptr<Shake> create_shake_backend() override;
    std::unique_ptr<NonbondedForce> create_nonbonded_force_backend() override;
};
