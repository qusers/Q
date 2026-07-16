#pragma once

#include <cstdint>
#include <memory>
#include <vector>

#include "host_device_buffer.h"
#include "shake.h"

// Reusable device-side solver for a set of independent molecules. Constraints
// for each molecule must be contiguous and in the desired serial order.
class CudaSerialConstraintSolver {
   public:
    void init(
        const std::vector<int>& molecule_constraint_counts,
        const std::vector<ShakeBond>& shake_bonds);
    void apply(coord_t* d_coords, const coord_t* d_xcoords, const double* d_winv);
    bool enabled() const { return n_constraints_ > 0; }

   private:
    int n_molecules_ = 0;
    int n_constraints_ = 0;
    std::unique_ptr<HostDeviceBuffer<int>> molecule_constraint_offsets_;
    std::unique_ptr<HostDeviceBuffer<ShakeBond>> shake_bonds_;
    std::unique_ptr<HostDeviceBuffer<std::uint8_t>> ready_;
    std::unique_ptr<HostDeviceBuffer<int>> failed_;
};

// GPU implementation of the original Q/CPU SHAKE update schedule.
//
// Molecules are independent, so they are processed in parallel.  Within one
// molecule a single CUDA thread visits constraints serially in topology order
// and retains the original persistent-ready behavior.
class CudaSerialShake final : public Shake {
   public:
    void apply(Context& ctx, HostDeviceBuffer<coord_t>& xcoords) override;
    void initial_shake(Context& ctx) override;

   protected:
    void init_backend(Context& ctx) override;

   private:
    void apply_to(Context& ctx, coord_t* d_coords, coord_t* d_xcoords);

    CudaSerialConstraintSolver solver_;
};
