#pragma once

#include "constraint_force.h"
#include "lincs.h"

struct CudaLincsBond {
    int ai;             // first atom, zero-based
    int aj;             // second atom, zero-based
    double target;      /// target bond length
    double mass_scale;  // 1 / sqrt(winv[ai] + winv[aj])
};

struct CudaLincsSmallData {  // For small groups
    int n_blocks = 0;
    int n_packed_threads = 0;
    int n_real_constraints = 0;
    int max_neighbors = 0;

    std::unique_ptr<HostDeviceBuffer<CudaLincsBond>> bonds;
    std::unique_ptr<HostDeviceBuffer<int>> neighbor_counts;

    std::unique_ptr<HostDeviceBuffer<int>> neighbor_local_indices;
    std::unique_ptr<HostDeviceBuffer<double>> mass_factors;

    std::unique_ptr<HostDeviceBuffer<double>> matrix_a;
    std::unique_ptr<HostDeviceBuffer<int>> failed_constraint;
};

class CudaLincs : public ConstraintForce {
   public:
    explicit CudaLincs(LincsSettings settings = {}) : lincs_settings_(settings) {};
    void apply(Context&, HostDeviceBuffer<coord_t>& xcoords) override;
    void initial_constraint(Context&) override;
    void cleanup() override;

   protected:
    void init_backend(Context&) override;

   private:
    void apply_to(Context& ctx, coord_t* coords, const coord_t* xcoords);
    LincsSettings lincs_settings_;
    CudaLincsSmallData lincs_small_data_;
};
