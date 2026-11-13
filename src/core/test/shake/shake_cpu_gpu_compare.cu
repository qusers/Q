#include <cuda_runtime.h>

#include <algorithm>
#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <memory>
#include <string>
#include <vector>

#include "common/include/constants.h"
#include "common/include/context.h"
#include "common/include/host_device_buffer.h"
#include "cpu/include/cpu_shake.h"
#include "cuda/include/cuda_shake_constraints.cuh"

namespace {

struct ShakeCase {
    std::string name;
    std::vector<coord_t> coords;
    std::vector<coord_t> xcoords;
    std::vector<real_t> winv;
    std::vector<int> mol_n_shakes;
    std::vector<shake_bond_t> shake_bonds;
};

template <typename T>
std::unique_ptr<HostDeviceBuffer<T>> make_buffer(const std::vector<T>& values) {
    auto buffer = std::make_unique<HostDeviceBuffer<T>>(values.size(), true, true);
    std::copy(values.begin(), values.end(), buffer->cpu_data_p);
    buffer->upload();
    return buffer;
}

bool constraints_conflict(const shake_bond_t& lhs, const shake_bond_t& rhs) {
    return lhs.ai == rhs.ai || lhs.ai == rhs.aj || lhs.aj == rhs.ai || lhs.aj == rhs.aj;
}

void build_shake_groups(
    const std::vector<shake_bond_t>& shake_bonds,
    const std::vector<int>& mol_n_shakes,
    std::vector<int>& mol_n_shake_groups,
    std::vector<int>& mol_shake_group_offset,
    std::vector<int>& shake_group_offset,
    std::vector<int>& shake_group_size,
    std::vector<int>& shake_group_indices) {
    int mol_shake_offset = 0;
    mol_n_shake_groups.assign(mol_n_shakes.size(), 0);
    mol_shake_group_offset.assign(mol_n_shakes.size(), 0);
    shake_group_offset.clear();
    shake_group_size.clear();
    shake_group_indices.clear();

    for (size_t mol = 0; mol < mol_n_shakes.size(); mol++) {
        std::vector<std::vector<int>> groups;
        for (int local_shake = 0; local_shake < mol_n_shakes[mol]; local_shake++) {
            const int shake_index = mol_shake_offset + local_shake;
            bool placed = false;
            for (auto& group : groups) {
                bool can_place = true;
                for (const int grouped_shake : group) {
                    if (constraints_conflict(shake_bonds[shake_index], shake_bonds[grouped_shake])) {
                        can_place = false;
                        break;
                    }
                }
                if (can_place) {
                    group.push_back(shake_index);
                    placed = true;
                    break;
                }
            }
            if (!placed) {
                groups.push_back({shake_index});
            }
        }

        mol_shake_group_offset[mol] = static_cast<int>(shake_group_offset.size());
        mol_n_shake_groups[mol] = static_cast<int>(groups.size());
        for (const auto& group : groups) {
            shake_group_offset.push_back(static_cast<int>(shake_group_indices.size()));
            shake_group_size.push_back(static_cast<int>(group.size()));
            shake_group_indices.insert(shake_group_indices.end(), group.begin(), group.end());
        }
        mol_shake_offset += mol_n_shakes[mol];
    }
}

void load_case(Context& ctx, const ShakeCase& test_case) {
    std::vector<int> mol_n_shake_groups;
    std::vector<int> mol_shake_group_offset;
    std::vector<int> shake_group_offset;
    std::vector<int> shake_group_size;
    std::vector<int> shake_group_indices;
    build_shake_groups(
        test_case.shake_bonds,
        test_case.mol_n_shakes,
        mol_n_shake_groups,
        mol_shake_group_offset,
        shake_group_offset,
        shake_group_size,
        shake_group_indices);

    ctx.command_info.requested_gpu = true;
    ctx.n_atoms = static_cast<int>(test_case.coords.size());
    ctx.n_molecules = static_cast<int>(test_case.mol_n_shakes.size());
    ctx.n_shake_constraints = static_cast<int>(test_case.shake_bonds.size());
    ctx.n_shake_groups = static_cast<int>(shake_group_offset.size());

    ctx.coords = make_buffer(test_case.coords);
    ctx.xcoords = make_buffer(test_case.xcoords);
    ctx.winv = make_buffer(test_case.winv);
    ctx.mol_n_shakes = make_buffer(test_case.mol_n_shakes);
    ctx.mol_n_shake_groups = make_buffer(mol_n_shake_groups);
    ctx.mol_shake_group_offset = make_buffer(mol_shake_group_offset);
    ctx.shake_group_offset = make_buffer(shake_group_offset);
    ctx.shake_group_size = make_buffer(shake_group_size);
    ctx.shake_group_indices = make_buffer(shake_group_indices);
    ctx.shake_bonds = make_buffer(test_case.shake_bonds);
}

std::vector<coord_t> run_cpu(Context& ctx, const ShakeCase& test_case) {
    load_case(ctx, test_case);
    calc_shake_constraints(ctx.coords->cpu_data_p, ctx.xcoords->cpu_data_p);
    return std::vector<coord_t>(ctx.coords->cpu_data_p, ctx.coords->cpu_data_p + ctx.n_atoms);
}

std::vector<coord_t> run_gpu(Context& ctx, const ShakeCase& test_case) {
    load_case(ctx, test_case);
    init_shake_constraints_kernel_data();
    calc_shake_constraints_host();
    cleanup_shake_constraints();
    return std::vector<coord_t>(ctx.coords->cpu_data_p, ctx.coords->cpu_data_p + ctx.n_atoms);
}

double max_abs_delta(const std::vector<coord_t>& lhs, const std::vector<coord_t>& rhs) {
    double max_delta = 0.0;
    for (size_t i = 0; i < lhs.size(); i++) {
        max_delta = std::max(max_delta, std::fabs(static_cast<double>(lhs[i].x) - static_cast<double>(rhs[i].x)));
        max_delta = std::max(max_delta, std::fabs(static_cast<double>(lhs[i].y) - static_cast<double>(rhs[i].y)));
        max_delta = std::max(max_delta, std::fabs(static_cast<double>(lhs[i].z) - static_cast<double>(rhs[i].z)));
    }
    return max_delta;
}

void print_coords(const char* label, const std::vector<coord_t>& coords) {
    std::printf("%s\n", label);
    for (size_t i = 0; i < coords.size(); i++) {
        std::printf(
            "  atom %zu: %.12g %.12g %.12g\n",
            i + 1,
            static_cast<double>(coords[i].x),
            static_cast<double>(coords[i].y),
            static_cast<double>(coords[i].z));
    }
}

ShakeCase make_iterative_case() {
    return {
        "iterative_shared_atom_constraints",
        {
            {static_cast<real_t>(0.0), static_cast<real_t>(0.0), static_cast<real_t>(0.0)},
            {static_cast<real_t>(1.16), static_cast<real_t>(0.08), static_cast<real_t>(0.0)},
            {static_cast<real_t>(2.05), static_cast<real_t>(-0.03), static_cast<real_t>(0.02)},
        },
        {
            {static_cast<real_t>(0.0), static_cast<real_t>(0.0), static_cast<real_t>(0.0)},
            {static_cast<real_t>(1.0), static_cast<real_t>(0.0), static_cast<real_t>(0.0)},
            {static_cast<real_t>(2.0), static_cast<real_t>(0.0), static_cast<real_t>(0.0)},
        },
        {static_cast<real_t>(1.0), static_cast<real_t>(1.0), static_cast<real_t>(1.0)},
        {2},
        {
            {1, 2, static_cast<real_t>(1.0), false},
            {2, 3, static_cast<real_t>(1.0), false},
        },
    };
}

ShakeCase make_ready_constraint_case() {
    const double current_dist = std::sqrt(1.0 - 0.5 * shake_tol);
    return {
        "ready_constraint_final_correction",
        {
            {static_cast<real_t>(0.0), static_cast<real_t>(0.0), static_cast<real_t>(0.0)},
            {static_cast<real_t>(current_dist), static_cast<real_t>(0.0), static_cast<real_t>(0.0)},
        },
        {
            {static_cast<real_t>(0.0), static_cast<real_t>(0.0), static_cast<real_t>(0.0)},
            {static_cast<real_t>(1.0), static_cast<real_t>(0.0), static_cast<real_t>(0.0)},
        },
        {static_cast<real_t>(1.0), static_cast<real_t>(1.0)},
        {1},
        {
            {1, 2, static_cast<real_t>(1.0), false},
        },
    };
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

    Context ctx;
    const std::vector<ShakeCase> cases = {
        make_iterative_case(),
        make_ready_constraint_case(),
    };

#ifdef QDYN_SPFP
    constexpr double coord_tolerance = 1.0e-6;
#else
    constexpr double coord_tolerance = 1.0e-12;
#endif

    bool failed = false;
    for (const auto& test_case : cases) {
        const std::vector<coord_t> cpu_coords = run_cpu(ctx, test_case);
        const std::vector<coord_t> gpu_coords = run_gpu(ctx, test_case);
        const double delta = max_abs_delta(cpu_coords, gpu_coords);

        if (delta > coord_tolerance) {
            failed = true;
            std::printf("[FAIL] %s: max_abs_delta %.12g > %.12g\n", test_case.name.c_str(), delta, coord_tolerance);
            print_coords("CPU SHAKE coords:", cpu_coords);
            print_coords("GPU SHAKE coords:", gpu_coords);
        } else {
            std::printf("[PASS] %s: max_abs_delta %.12g\n", test_case.name.c_str(), delta);
        }
    }

    return failed ? EXIT_FAILURE : EXIT_SUCCESS;
}
