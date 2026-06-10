#pragma once

#include "handler.h"


class CpuHandler : public Handler {
   public:
    static CpuHandler& instance() {
        static CpuHandler h;
        return h;
    }

    void shutdown() override;

   protected:
    void initialize_backend() override;
    void calc_internal_forces(int iteration) override;
    void calc_nonbonded_forces() override;
    void calc_temperature() override;
    void calc_leapfrog() override;
    CpuHandler() = default;

    std::unique_ptr<Shake> create_shake_backend() override;
};
