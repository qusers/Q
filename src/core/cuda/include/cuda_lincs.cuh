#pragma once

#include "constraint_force.h"
#include "lincs.h"

struct CudaLincsBond {
    int ai;             // first atom, zero-based
    int aj;             // second atom, zero-based
    double target;      /// target bond length
    double mass_scale;  // 1 / sqrt(winv[ai] + winv[aj])
};

struct CudaLincsData {
    int n_blocks = 0;
    int n_constraints = 0;
    int max_neighbors = 0;
    int max_work_unit_size = 0;

    std::unique_ptr<HostDeviceBuffer<int>> work_unit_offsets;
    std::unique_ptr<HostDeviceBuffer<CudaLincsBond>> bonds;
    std::unique_ptr<HostDeviceBuffer<int>> neighbor_counts;
    std::unique_ptr<HostDeviceBuffer<int>> neighbor_indices;
    std::unique_ptr<HostDeviceBuffer<double>> mass_factors;

    std::unique_ptr<HostDeviceBuffer<coord_t>> directions;
    std::unique_ptr<HostDeviceBuffer<double>> matrix_a;
    std::unique_ptr<HostDeviceBuffer<double>> q0;
    std::unique_ptr<HostDeviceBuffer<double>> q1;
    std::unique_ptr<HostDeviceBuffer<double>> solution;

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
    CudaLincsData lincs_data_;
};
