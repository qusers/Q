#pragma once

#include "context.h"

class Handler {
   public:
    virtual ~Handler() = default;
    virtual void initialize() = 0;

    virtual void shutdown() = 0;

    void run_iteration(int iteration);

    void run(int num_iterations);

   protected:
    Handler() = default;
    Handler(const Handler&) = delete;
    Handler& operator=(const Handler&) = delete;

    bool initialized_ = false;
    Context& ctx = Context::instance();

    virtual void calc_internal_forces(int iteration) = 0;
    virtual void calc_nonbonded_forces() = 0;
    virtual void calc_temperature() = 0;
    virtual void calc_leapfrog() = 0;


    void update_energy_totals();
    void print_outputs(int iteration);
    virtual void reset_energies();
};
