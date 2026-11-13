#pragma once
#include <array>
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
    void initialize(const CommandInfo& command);        // CLI: parse .inp, system + first stage
    void initialize_system(const ParseResult& parsed);  // once. System fields, plus the starting coords/velocities from the restart file if there is one
    void configure(const StageSettings& stage);         // per stage
    void validate_stage(const StageSettings& stage) const { ctx.validate_stage(stage); }  // throws, changes nothing
    void set_state(const SystemState& s);
    SystemState get_state();
    void run();  // runs the configured stage

   protected:
    explicit Handler(bool use_gpu) { ctx.command_info.requested_gpu = use_gpu; }
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
    void start_fresh(int seed, double initial_temperature);
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
    void prepare_offdiags(int iteration);
    std::unique_ptr<ConstraintForce> constraint_force_;
    std::unique_ptr<NonbondedForce> nonbonded_force_;
    std::unique_ptr<BondedForce> bonded_force_;
    std::unique_ptr<RestraintForce> restraint_force_;
    std::unique_ptr<Temperature> temperature_;
    std::unique_ptr<Integrator> integrator_;
    std::unique_ptr<WaterBoundaryForce> water_boundary_force_;
    std::vector<std::unique_ptr<BaseOutput>> outputs_;
    Context ctx;

    std::vector<double> current_theta_corr();
    void rebuild_constraints_and_bonded();

    ParseResult system_;
    bool system_loaded_ = false;
    bool configured_ = false;
    bool state_ready_ = false;
    bool broken_ = false;                     // NEW: set when run() throws
    std::array<bool, 3> constraint_flags_{};  // solvent, solute, hydrogens of the last build
};
