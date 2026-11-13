// Parity test for CudaSerialConstraintSolver.
//
// The serial CUDA solver promises to reproduce the original Q SHAKE update
// schedule (src/q6/md.f90 shake(), mirrored by CpuShake::apply_to):
//   - constraints of one molecule are visited serially in topology order;
//   - a constraint whose relative error drops below shake_tol is marked ready
//     and receives that visit's correction, then is skipped for the rest of
//     the molecule's iteration loop, even if later corrections move it again;
//   - a molecule that is not ready after shake_max_iter sweeps is a failure.
//
// The test runs the solver on small systems and compares the coordinates with
// a host implementation of that schedule. It does not require every bond to be
// within tolerance afterwards: with persistent ready flags Q itself does not
// guarantee that, and the "coupled_ready_skip" case checks that this
// behavior is kept.

#include <cuda_runtime.h>

#include <algorithm>
#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <string>
#include <vector>

#include "common/include/constants.h"
#include "common/include/host_device_buffer.h"
#include "common/include/qgpu_error.h"
#include "cuda/include/cuda_serial_shake.cuh"

namespace {

// Both sides do the same double-precision operations in the same order
// (device code is built with --fmad=false), so they should agree to roundoff.
constexpr double kCoordTolerance = 1.0e-12;

struct ShakeCase {
    std::string name;
    std::vector<coord_t> coords;
    std::vector<coord_t> xcoords;
    std::vector<double> winv;
    std::vector<int> molecule_constraint_counts;
    std::vector<ConstraintBond> bonds;
};

// Host version of the Q SHAKE loop. Returns false if a molecule does not
// converge within shake_max_iter sweeps.
bool reference_shake(const ShakeCase& c, std::vector<coord_t>& x) {
    const std::vector<coord_t>& xx = c.xcoords;
    bool all_converged = true;
    int first = 0;
    for (const int count : c.molecule_constraint_counts) {
        const int end = first + count;
        std::vector<bool> ready(count, false);
        bool converged = false;
        int nits = 0;
        do {
            converged = true;
            for (int ic = first; ic < end; ic++) {
                if (ready[ic - first]) continue;
                const ConstraintBond& b = c.bonds[ic];
                const int i = b.ai - 1;
                const int j = b.aj - 1;
                const double xij_x = x[i].x - x[j].x;
                const double xij_y = x[i].y - x[j].y;
                const double xij_z = x[i].z - x[j].z;
                const double xij2 = xij_x * xij_x + xij_y * xij_y + xij_z * xij_z;
                const double diff = b.dist2 - xij2;
                if (std::fabs(diff) < shake_tol * b.dist2) {
                    ready[ic - first] = true;
                } else {
                    converged = false;
                }
                const double xxij_x = xx[i].x - xx[j].x;
                const double xxij_y = xx[i].y - xx[j].y;
                const double xxij_z = xx[i].z - xx[j].z;
                const double scp = xij_x * xxij_x + xij_y * xxij_y + xij_z * xxij_z;
                const double corr = diff / (2.0 * scp * (c.winv[i] + c.winv[j]));
                x[i].x += xxij_x * corr * c.winv[i];
                x[i].y += xxij_y * corr * c.winv[i];
                x[i].z += xxij_z * corr * c.winv[i];
                x[j].x -= xxij_x * corr * c.winv[j];
                x[j].y -= xxij_y * corr * c.winv[j];
                x[j].z -= xxij_z * corr * c.winv[j];
            }
            nits++;
        } while (nits < shake_max_iter && !converged);
        if (!converged) all_converged = false;
        first = end;
    }
    return all_converged;
}

struct GpuResult {
    std::vector<coord_t> coords;
    bool threw_constraint_error = false;
};

GpuResult run_gpu(const ShakeCase& c) {
    auto coords = HostDeviceBuffer<coord_t>::from_vector(c.coords, true);
    auto xcoords = HostDeviceBuffer<coord_t>::from_vector(c.xcoords, true);
    auto winv = HostDeviceBuffer<double>::from_vector(c.winv, true);

    GpuResult result;
    CudaSerialConstraintSolver solver;
    solver.init(c.molecule_constraint_counts, c.bonds);
    try {
        solver.apply(coords->gpu_data_p, xcoords->gpu_data_p, winv->gpu_data_p);
    } catch (const ConstraintError&) {
        result.threw_constraint_error = true;
    }
    coords->download();
    result.coords.assign(coords->cpu_data_p, coords->cpu_data_p + c.coords.size());
    return result;
}

double max_abs_delta(const std::vector<coord_t>& a, const std::vector<coord_t>& b) {
    double delta = 0.0;
    for (size_t i = 0; i < a.size(); i++) {
        delta = std::max(delta, std::fabs(a[i].x - b[i].x));
        delta = std::max(delta, std::fabs(a[i].y - b[i].y));
        delta = std::max(delta, std::fabs(a[i].z - b[i].z));
    }
    return delta;
}

double relative_error(const std::vector<coord_t>& x, const ConstraintBond& b) {
    const coord_t& p = x[b.ai - 1];
    const coord_t& q = x[b.aj - 1];
    const double dx = p.x - q.x;
    const double dy = p.y - q.y;
    const double dz = p.z - q.z;
    return std::fabs(dx * dx + dy * dy + dz * dz - b.dist2) / b.dist2;
}

// Two constraints share atom 2. Bond 1-2 starts exactly at its target, so it is
// marked ready on the first visit; correcting 2-3 then moves atom 2, and Q
// never revisits 1-2. Q leaves 1-2 with a relative error near 1e-1.
ShakeCase coupled_ready_skip_case() {
    return {
        "coupled_ready_skip",
        {{-0.3, 0.0, 0.0}, {0.7, 0.0, 0.0}, {1.6, 0.0, 0.0}},
        {{0.0, 0.0, 0.0}, {1.0, 0.0, 0.0}, {2.0, 0.0, 0.0}},
        {1.0, 1.0, 1.0},
        {2},
        {{1, 2, 1.0}, {2, 3, 1.0}},
    };
}

// Two shared-atom constraints that both start out of tolerance and need
// several sweeps.
ShakeCase iterative_shared_atom_case() {
    return {
        "iterative_shared_atom",
        {{0.0, 0.0, 0.0}, {1.16, 0.08, 0.0}, {2.05, -0.03, 0.02}},
        {{0.0, 0.0, 0.0}, {1.0, 0.0, 0.0}, {2.0, 0.0, 0.0}},
        {1.0, 0.5, 2.0},
        {2},
        {{1, 2, 1.0}, {2, 3, 1.0}},
    };
}

// A constraint already within tolerance still receives the correction of the
// visit on which it becomes ready.
ShakeCase ready_final_correction_case() {
    const double current_dist = std::sqrt(1.0 - 0.5 * shake_tol);
    return {
        "ready_final_correction",
        {{0.0, 0.0, 0.0}, {current_dist, 0.0, 0.0}},
        {{0.0, 0.0, 0.0}, {1.0, 0.0, 0.0}},
        {1.0, 1.0},
        {1},
        {{1, 2, 1.0}},
    };
}

// The three cases above as separate molecules of one system, with an
// unconstrained molecule in between, so one kernel launch handles several
// molecules and the per-molecule offsets.
ShakeCase multi_molecule_case() {
    const ShakeCase parts[] = {
        coupled_ready_skip_case(),
        iterative_shared_atom_case(),
        ready_final_correction_case(),
    };
    ShakeCase c;
    c.name = "multi_molecule";
    for (int p = 0; p < 3; p++) {
        const int atom_offset = static_cast<int>(c.coords.size());
        c.coords.insert(c.coords.end(), parts[p].coords.begin(), parts[p].coords.end());
        c.xcoords.insert(c.xcoords.end(), parts[p].xcoords.begin(), parts[p].xcoords.end());
        c.winv.insert(c.winv.end(), parts[p].winv.begin(), parts[p].winv.end());
        c.molecule_constraint_counts.push_back(parts[p].molecule_constraint_counts[0]);
        for (ConstraintBond b : parts[p].bonds) {
            b.ai += atom_offset;
            b.aj += atom_offset;
            c.bonds.push_back(b);
        }
        if (p == 0) {
            c.coords.push_back({5.0, 5.0, 5.0});
            c.xcoords.push_back({5.0, 5.0, 5.0});
            c.winv.push_back(1.0);
            c.molecule_constraint_counts.push_back(0);
        }
    }
    return c;
}

// The old and new bond vectors are perpendicular, so the correction is
// infinite and the molecule cannot converge. The solver must report this with
// ConstraintError.
ShakeCase non_converging_case() {
    return {
        "non_converging",
        {{0.0, 0.0, 0.0}, {1.1, 0.0, 0.0}},
        {{0.0, 0.0, 0.0}, {0.0, 1.0, 0.0}},
        {1.0, 1.0},
        {1},
        {{1, 2, 1.0}},
    };
}

bool check_parity(const ShakeCase& c, std::vector<coord_t>* reference_out = nullptr) {
    std::vector<coord_t> reference = c.coords;
    if (!reference_shake(c, reference)) {
        std::printf("[FAIL] %s: host reference did not converge\n", c.name.c_str());
        return false;
    }
    const GpuResult gpu = run_gpu(c);
    if (gpu.threw_constraint_error) {
        std::printf("[FAIL] %s: solver threw ConstraintError\n", c.name.c_str());
        return false;
    }
    const double delta = max_abs_delta(gpu.coords, reference);
    if (!std::isfinite(delta) || delta > kCoordTolerance) {
        std::printf(
            "[FAIL] %s: max |GPU - Q reference| %.12g > %.12g\n",
            c.name.c_str(), delta, kCoordTolerance);
        return false;
    }
    std::printf("[PASS] %s: max |GPU - Q reference| %.12g\n", c.name.c_str(), delta);
    if (reference_out) *reference_out = reference;
    return true;
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

    bool ok = true;

    // Parity alone could pass if the reference were wrong in the same way as
    // the kernel, so the first and third cases also check that the reference
    // really shows the persistent-ready and final-correction behavior.
    {
        const ShakeCase c = coupled_ready_skip_case();
        std::vector<coord_t> reference;
        if (check_parity(c, &reference)) {
            const double err = relative_error(reference, c.bonds[0]);
            if (err < shake_tol) {
                std::printf(
                    "[FAIL] %s: bond 1-2 relative error %.12g; expected it to stay above %g\n",
                    c.name.c_str(), err, static_cast<double>(shake_tol));
                ok = false;
            }
        } else {
            ok = false;
        }
    }

    ok = check_parity(iterative_shared_atom_case()) && ok;

    {
        const ShakeCase c = ready_final_correction_case();
        std::vector<coord_t> reference;
        if (check_parity(c, &reference)) {
            if (max_abs_delta(reference, c.coords) == 0.0) {
                std::printf("[FAIL] %s: expected the final correction to move the atoms\n", c.name.c_str());
                ok = false;
            }
        } else {
            ok = false;
        }
    }

    ok = check_parity(multi_molecule_case()) && ok;

    {
        const ShakeCase c = non_converging_case();
        std::vector<coord_t> reference = c.coords;
        const bool reference_converged = reference_shake(c, reference);
        const GpuResult gpu = run_gpu(c);
        if (reference_converged || !gpu.threw_constraint_error) {
            std::printf(
                "[FAIL] %s: reference converged=%d, solver threw ConstraintError=%d\n",
                c.name.c_str(), reference_converged, gpu.threw_constraint_error);
            ok = false;
        } else {
            std::printf("[PASS] %s: solver threw ConstraintError\n", c.name.c_str());
        }
    }

    return ok ? EXIT_SUCCESS : EXIT_FAILURE;
}
