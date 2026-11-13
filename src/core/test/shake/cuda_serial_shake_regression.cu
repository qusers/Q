#include <cuda_runtime.h>

#include <algorithm>
#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <vector>

#include "common/include/constants.h"
#include "common/include/host_device_buffer.h"
#include "cuda/include/cuda_serial_shake.cuh"

namespace {

double relative_constraint_error(
    const coord_t* coords,
    const ConstraintBond& bond) {
    const int ai = bond.ai - 1;
    const int aj = bond.aj - 1;
    const double dx = coords[ai].x - coords[aj].x;
    const double dy = coords[ai].y - coords[aj].y;
    const double dz = coords[ai].z - coords[aj].z;
    const double dist2 = dx * dx + dy * dy + dz * dz;
    return std::fabs(dist2 - bond.dist2) / bond.dist2;
}

int skip_without_cuda_device() {
    int device_count = 0;
    const cudaError_t status = cudaGetDeviceCount(&device_count);
    if (status != cudaSuccess || device_count == 0) {
        std::printf("SKIP: CUDA device is not available: %s\n", cudaGetErrorString(status));
        return 77;
    }
    return EXIT_SUCCESS;
}

}  // namespace

int main() {
    const int cuda_status = skip_without_cuda_device();
    if (cuda_status != EXIT_SUCCESS) return cuda_status;

    // Two coupled constraints share atom 2. The first bond starts exactly at
    // its target and the second does not. Persistent-ready SHAKE permanently
    // skipped the first bond after correcting the second and returned with its
    // relative squared-distance error near 1e-1.
    const std::vector<coord_t> initial_coords = {
        {-0.3, 0.0, 0.0},
        {0.7, 0.0, 0.0},
        {1.6, 0.0, 0.0},
    };
    const std::vector<coord_t> reference_coords = {
        {0.0, 0.0, 0.0},
        {1.0, 0.0, 0.0},
        {2.0, 0.0, 0.0},
    };
    const std::vector<double> inverse_masses = {1.0, 1.0, 1.0};
    const std::vector<ConstraintBond> constraints = {
        {1, 2, 1.0},
        {2, 3, 1.0},
    };

    auto coords = HostDeviceBuffer<coord_t>::from_vector(initial_coords, true);
    auto xcoords = HostDeviceBuffer<coord_t>::from_vector(reference_coords, true);
    auto winv = HostDeviceBuffer<double>::from_vector(inverse_masses, true);

    CudaSerialConstraintSolver solver;
    solver.init({static_cast<int>(constraints.size())}, constraints);
    solver.apply(coords->gpu_data_p, xcoords->gpu_data_p, winv->gpu_data_p);
    coords->download();

    double max_error = 0.0;
    for (const ConstraintBond& constraint : constraints) {
        max_error = std::max(
            max_error,
            relative_constraint_error(coords->cpu_data_p, constraint));
    }

    if (!std::isfinite(max_error) || max_error >= shake_tol) {
        std::printf(
            "[FAIL] coupled constraints: max relative error %.12g >= %.12g\n",
            max_error,
            static_cast<double>(shake_tol));
        return EXIT_FAILURE;
    }

    std::printf(
        "[PASS] coupled constraints: max relative error %.12g\n",
        max_error);
    return EXIT_SUCCESS;
}
