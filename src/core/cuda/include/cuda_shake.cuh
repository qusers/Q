#pragma once
#include "constraint_force.h"
#include "cuda_serial_shake.cuh"
#include "cuda_settle.cuh"
#include "host_device_buffer.h"

struct ShakeNetwork {
    int center;
    int n_hydrogens;
    int hydrogens[3];  // at most 3
    double dist2[3];
    double center_winv;
    double hydrogen_winv[3];
};

class CudaShake final : public ConstraintForce {
   public:
    explicit CudaShake(bool serial_q_molecules = false)
        : serial_q_molecules_(serial_q_molecules) {}

    void apply(Context& ctx, HostDeviceBuffer<coord_t>& xcoords) override;
    void initial_constraint(Context& ctx) override;
    void cleanup() override;
    void init_from_bonds(Context& ctx, const std::vector<ConstraintBond>& bonds) override;

   protected:
    void init_backend(Context& ctx) override;

   private:
    void apply_to(Context& ctx, coord_t* d_coords, coord_t* d_xcoords);

    bool is_init_backend = false;
    bool serial_q_molecules_ = false;
    CudaSerialConstraintSolver serial_q_solver_;
    // CudaSettleFastWaterSolver settle_fast_water_solver_;

    std::unique_ptr<HostDeviceBuffer<ShakeNetwork>> shake_networks;

    std::unique_ptr<HostDeviceBuffer<ConstraintBond>> fallback_shake_bonds;
    std::unique_ptr<HostDeviceBuffer<int>> fallback_color_offsets;
    std::unique_ptr<HostDeviceBuffer<int>> fallback_unconverged;
    std::unique_ptr<HostDeviceBuffer<int>> shake_network_failed;
    int fallback_n_colors = 0;
    int fallback_coop_blocks = 0;

    void find_shake_network(Context& ctx, const std::vector<ConstraintBond>& bonds, std::vector<bool>& optimized);
    void find_fallback_shake_bond(Context& ctx, const std::vector<ConstraintBond>& bonds, std::vector<bool>& optimized);
    void find_serial_q_molecule_bonds(Context& ctx, const std::vector<ConstraintBond>& bonds, std::vector<bool>& optimized);
};
