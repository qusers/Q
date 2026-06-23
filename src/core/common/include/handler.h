#pragma once

#include <memory>
#include <vector>

#include "base_output.h"
#include "context.h"
#include "nonbonded_force.h"
#include "shake.h"

class Handler {
   public:
    Shake& shake();
    NonbondedForce& nonbonded_force();
    virtual ~Handler() = default;
    virtual void initialize();

    virtual void shutdown() = 0;

    void run_iteration(int iteration);

    void run(int num_iterations);

   protected:
    Handler() = default;
    Handler(const Handler&) = delete;
    Handler& operator=(const Handler&) = delete;

    bool initialized_ = false;
    Context& ctx = Context::instance();
    std::vector<std::unique_ptr<BaseOutput>> outputs_;

    virtual void calc_internal_forces(int iteration) = 0;
    virtual void calc_nonbonded_forces() = 0;
    virtual void calc_temperature() = 0;
    virtual void calc_leapfrog() = 0;
    virtual void initialize_backend() = 0;

    void calc_final_potential(int iteration);

    void update_energy_totals();
    void print_outputs(int iteration);
    void create_outputs();
    void init_outputs();
    void finish_outputs();
    void shutdown_outputs();
    virtual void reset_energies();

    std::unique_ptr<Shake> shake_;
    std::unique_ptr<NonbondedForce> nonbonded_force_;
    virtual std::unique_ptr<Shake> create_shake_backend() = 0;
    virtual std::unique_ptr<NonbondedForce> create_nonbonded_force_backend() = 0;
};
