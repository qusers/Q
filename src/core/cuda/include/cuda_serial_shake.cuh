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
        const std::vector<ShakeBond>& shake_bonds, int n_atoms, int max_replicates);
    void apply(coord_t* d_coords, const coord_t* d_xcoords, const double* d_winv, int n_replicates);
    bool enabled() const { return n_constraints_ > 0; }

   private:
    int n_molecules_ = 0;
    int n_constraints_ = 0;
    int n_atoms_ = 0;
    int max_replicates_ = 0;
    std::unique_ptr<HostDeviceBuffer<int>> molecule_constraint_offsets_;
    std::unique_ptr<HostDeviceBuffer<ShakeBond>> shake_bonds_;
    std::unique_ptr<HostDeviceBuffer<std::uint8_t>> ready_;
    std::unique_ptr<HostDeviceBuffer<int>> failed_;
};
