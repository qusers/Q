#pragma once

#include <memory>
#include <string>
#include <vector>

#include "base_output.h"
#include "context.h"

class Handler {
   public:
    virtual ~Handler() = default;
    virtual void initialize();

    virtual void shutdown() = 0;

    void run_iteration(int iteration);

    void run(int num_iterations);
    void run_final_iteration(int iteration);

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
    virtual void prepare_force_dump();


    void update_energy_totals();
    void dump_force_group(const std::string& group, int iteration);
    void print_outputs(int iteration);
    void print_final_outputs(int iteration);
    void create_outputs();
    void init_outputs();
    void finish_outputs();
    void shutdown_outputs();
    virtual void reset_energies();
};
