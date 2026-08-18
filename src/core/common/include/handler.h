#pragma once

#include <memory>
#include <vector>

#include "base_output.h"
#include "bonded_force.h"
#include "constraint_force.h"
#include "context.h"
#include "integrator.h"
#include "nonbonded_force.h"
#include "restraint_force.h"
#include "temperature.h"
#include "water_boundary_force.h"

class Handler {
   public:
    virtual ~Handler() = default;
    void initialize(const CommandInfo& command);
    void run();

   protected:
    Handler() = default;
    Handler(const Handler&) = delete;
    Handler& operator=(const Handler&) = delete;

    virtual std::unique_ptr<ConstraintForce> create_constraint_force_backend() = 0;
    virtual std::unique_ptr<NonbondedForce> create_nonbonded_force_backend() = 0;
    virtual std::unique_ptr<BondedForce> create_bonded_force_backend() = 0;
    virtual std::unique_ptr<RestraintForce> create_restraint_force_backend() = 0;
    virtual std::unique_ptr<Temperature> create_temperature_backend() = 0;
    virtual std::unique_ptr<Integrator> create_integrator_backend() = 0;
    virtual std::unique_ptr<WaterBoundaryForce> create_water_boundary_force_backend() = 0;

   private:
    void stop_cm_translation();
    void calc_internal_forces(int iteration);
    void calc_nonbonded_forces();
    void calc_temperature();
    void calc_leapfrog();

    void calc_final_potential(int iteration);
    void update_energy_totals();
    void print_outputs(int iteration);
    void create_outputs();
    void init_outputs();
    void finish_outputs();
    void shutdown_outputs();
    void reset_energies();
    void run_iteration(int iteration);
    std::unique_ptr<ConstraintForce> constraint_force_;
    std::unique_ptr<NonbondedForce> nonbonded_force_;
    std::unique_ptr<BondedForce> bonded_force_;
    std::unique_ptr<RestraintForce> restraint_force_;
    std::unique_ptr<Temperature> temperature_;
    std::unique_ptr<Integrator> integrator_;
    std::unique_ptr<WaterBoundaryForce> water_boundary_force_;
    std::vector<std::unique_ptr<BaseOutput>> outputs_;
    Context ctx;
};
