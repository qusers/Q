#include <cooperative_groups.h>

#include <queue>
#include <stdexcept>
#include <vector>

#include "cuda_lincs.cuh"
#include "geometry.h"
#include "union_find.h"

namespace cg = cooperative_groups;

namespace {
const int BLOCK_THREAD = 256;

int get_same_atom_in_two_bond(const ConstraintBond& b1, const ConstraintBond& b2) {
    if (b1.ai == b2.ai || b1.ai == b2.aj) {
        return b1.ai;
    } else {
        return b1.aj;
    }
}

int get_sign(const int a, const ConstraintBond& bond) {
    // we define the direction is bond.ai - bond.aj
    if (a == bond.ai) return 1;
    return -1;
}

__device__ void solve(
    int n_packed_threads,  // total threads
    int packed,            // current thread idx in the whole kernel
    int local,             // current local thread idx
    int neighbor_count,    // neighbor count
    bool is_dummy,
    const CudaLincsBond& bond,
    const coord_t& direction,

    int expansion_order,

    const int* neighbor_local_indices,
    double* matrix_a,
    double* shared_rhs,

    const double* winv,
    coord_t* coords,
    const coord_t* reference,
    int* failed) {
    double solution = shared_rhs[local];
    /*
     * Neumann expansion:
     *
     * solution =
     *     rhs + A*rhs + A^2*rhs + ... + A^order*rhs
     */

    for (int i = 0; i < expansion_order; i++) {
        __syncthreads();

        const int input_buffer = i & 1;
        const int output_buffer = input_buffer ^ 1;

        double next_value = 0;

        for (int j = 0; j < neighbor_count; j++) {
            const int storage_idx = j * n_packed_threads + packed;
            const int neighbor_local = neighbor_local_indices[storage_idx];
            next_value += matrix_a[storage_idx] * shared_rhs[input_buffer * BLOCK_THREAD + neighbor_local];
        }

        shared_rhs[output_buffer * BLOCK_THREAD + local] = next_value;
        /*
        If it's for big group, we should wait until all block finishs.

        */
        solution += next_value;
    }

    /*
     * Apply the initial coordinate correction:
     *
     * lambda = mass_scale * solution
     *
     * ai -= winv[ai] * lambda * direction
     * aj += winv[aj] * lambda * direction
     */
    if (!is_dummy) {
        const double lambda = bond.mass_scale * solution;
        if (isfinite(lambda)) {
            atomicAdd(&coords[bond.ai].x, -winv[bond.ai] * lambda * direction.x);
            atomicAdd(&coords[bond.ai].y, -winv[bond.ai] * lambda * direction.y);
            atomicAdd(&coords[bond.ai].z, -winv[bond.ai] * lambda * direction.z);
            atomicAdd(&coords[bond.aj].x, winv[bond.aj] * lambda * direction.x);
            atomicAdd(&coords[bond.aj].y, winv[bond.aj] * lambda * direction.y);
            atomicAdd(&coords[bond.aj].z, winv[bond.aj] * lambda * direction.z);
        } else {
            atomicCAS(failed, -1, packed);
        }
    }
}

__global__ void small_lincs_kernel(
    const int n_packed_threads,
    const int expansion_order,
    const int minimum_rotation_iterations,
    const int maximum_rotation_iterations,
    const double accuracy_tolerance,
    const CudaLincsBond* bonds,
    const int* neighbor_counts,
    const int* neighbor_local_indices,
    const double* mass_factors,
    double* matrix_a,
    const double* winv,
    coord_t* coords,
    const coord_t* reference,
    int* failed

) {
    const int packed = blockIdx.x * blockDim.x + threadIdx.x;
    const int local = threadIdx.x;

    if (packed >= n_packed_threads) return;

    /*
     * Shared memory layout:
     *
     * directions[local]:
     *     old/reference unit direction of this constraint.
     *
     * rhs[0 * BLOCK_THREAD + local]:
     * rhs[1 * BLOCK_THREAD + local]:
     *     ping-pong vectors for the Neumann expansion
     */
    __shared__ coord_t shared_directions[BLOCK_THREAD];
    __shared__ double shared_rhs[2 * BLOCK_THREAD];
    __shared__ int block_ready;

    const CudaLincsBond& bond = bonds[packed];
    const bool is_dummy = bond.ai < 0;
    coord_t direction = {0, 0, 0};

    if (!is_dummy) {
        const coord_t old_bond = reference[bond.ai] - reference[bond.aj];
        const double old_d2 = norm2(old_bond);

        if (isfinite(old_d2) && old_d2 > 0) {
            const double inverse_length = rsqrt(old_d2);
            direction = old_bond * inverse_length;
        } else {
            atomicCAS(failed, -1, packed);
        }
    }

    shared_directions[local] = direction;
    __syncthreads();

    // Construct matrix-A
    const int count = neighbor_counts[packed];
    for (int i = 0; i < count; i++) {
        int storage_idx = i * n_packed_threads + packed;
        const int neighbor_local = neighbor_local_indices[storage_idx];
        matrix_a[storage_idx] = mass_factors[storage_idx] * dot(direction, shared_directions[neighbor_local]);
    }

    double rhs = 0;

    if (!is_dummy) {
        const coord_t current_bond = coords[bond.ai] - coords[bond.aj];
        rhs = bond.mass_scale * (dot(direction, current_bond) - bond.target);

        if (!isfinite(rhs)) {
            atomicCAS(failed, -1, packed);
            rhs = 0;
        }
    }
    shared_rhs[local] = rhs;
    solve(n_packed_threads, packed, local, count, is_dummy, bond, direction, expansion_order, neighbor_local_indices, matrix_a, shared_rhs, winv, coords, reference, failed);

    /*
     * Rotational/centripetal corrections.
     */
    for (int iteration = 0; iteration < maximum_rotation_iterations; iteration++) {
        __syncthreads();
        double rotation_rhs = 0;

        if (!is_dummy) {
            const coord_t current_bond = coords[bond.ai] - coords[bond.aj];
            const double current_length2 = norm2(current_bond);

            const double target_length2 = bond.target * bond.target;

            const double radicand = 2.0 * target_length2 - current_length2;

            if (isfinite(radicand) && radicand > 0) {
                rotation_rhs = bond.mass_scale * (bond.target - sqrt(radicand));
                if (!isfinite(rotation_rhs)) {
                    rotation_rhs = 0;
                    atomicCAS(failed, -1, packed);
                }
            } else {
                rotation_rhs = 0;
                atomicCAS(failed, -1, packed);
            }
        }
        shared_rhs[local] = rotation_rhs;
        solve(n_packed_threads, packed, local, count, is_dummy, bond, direction, expansion_order, neighbor_local_indices, matrix_a, shared_rhs, winv, coords, reference, failed);

        __syncthreads();
        if (local == 0) {
            block_ready = 1;
        }
        __syncthreads();
        if (!is_dummy) {
            const coord_t corrected_bond = coords[bond.ai] - coords[bond.aj];
            const double corrected_length = norm(corrected_bond);

            const double relative_error = fabs(corrected_length / bond.target - 1.0);
            if (!isfinite(corrected_length) || !isfinite(relative_error) || relative_error > accuracy_tolerance) {
                atomicExch(&block_ready, 0);
            }
        }
        __syncthreads();
        if (iteration + 1 >= minimum_rotation_iterations && block_ready != 0) {
            break;
        }
    }

    /*
     * Final validation. The preceding convergence check ended with
     * a block barrier, so all coordinate corrections are visible.
     */

    __syncthreads();
    if (!is_dummy) {
        const coord_t corrected_bond = coords[bond.ai] - coords[bond.aj];
        const double corrected_length = norm(corrected_bond);

        const double relative_error = fabs(corrected_length / bond.target - 1.0);
        if (!isfinite(corrected_length) || !isfinite(relative_error) || relative_error > accuracy_tolerance) {
            printf(">>> Lincs failed, i = %d, j = %d, d = %f, d0 = %f\n", bond.ai, bond.aj, corrected_length, bond.target);
            atomicCAS(failed, -1, packed);
        }
    }
}
// __restrict__: This pointer does not alias/overlap with the memory accessed through the other restricted pointers.
__global__ void big_lincs_prepare_direction_rhs_kernel(
    const int n_real_constraints,
    const CudaLincsBond* __restrict__ bonds,  //
    coord_t* __restrict__ directions,
    double* __restrict__ rhs_ping,
    double* __restrict__ solution,
    const coord_t* __restrict__ coords,
    const coord_t* __restrict__ reference,
    int* failed) {
    const int packed = blockIdx.x * blockDim.x + threadIdx.x;

    if (packed >= n_real_constraints) return;

    const CudaLincsBond& bond = bonds[packed];

    const coord_t old_bond = reference[bond.ai] - reference[bond.aj];

    const double old_length2 = norm2(old_bond);

    coord_t direction = {0, 0, 0};
    bool direction_valid = isfinite(old_length2) && old_length2 > 0;

    if (direction_valid) {
        direction = old_bond * rsqrt(old_length2);
    } else {
        atomicCAS(failed, -1, packed);
    }

    directions[packed] = direction;

    double rhs = 0;

    if (direction_valid) {
        const coord_t current_bond = coords[bond.ai] - coords[bond.aj];

        rhs = bond.mass_scale * (dot(direction, current_bond) - bond.target);

        if (!isfinite(rhs)) {
            rhs = 0;
            atomicCAS(failed, -1, packed);
        }
    }

    rhs_ping[packed] = rhs;
    solution[packed] = rhs;
}

__global__ void big_lincs_prepare_matrix_kernel(
    const int n_packed_threads,
    const int n_real_constraints,
    const int* __restrict__ neighbor_counts,
    const int* __restrict__ neighbor_indices,
    const double* __restrict__ mass_factors,
    const coord_t* __restrict__ directions,
    double* __restrict__ matrix_a,
    int* failed) {
    const int packed = blockIdx.x * blockDim.x + threadIdx.x;

    if (packed >= n_real_constraints) return;

    const coord_t direction = directions[packed];
    const int neighbor_count = neighbor_counts[packed];

    for (int slot = 0; slot < neighbor_count; slot++) {
        const int storage_idx = slot * n_packed_threads + packed;

        const int neighbor_idx = neighbor_indices[storage_idx];

        const double value = mass_factors[storage_idx] * dot(direction, directions[neighbor_idx]);

        if (isfinite(value)) {
            matrix_a[storage_idx] = value;
        } else {
            matrix_a[storage_idx] = 0;
            atomicCAS(failed, -1, packed);
        }
    }
}

__device__ void big_lincs_solve_grid(
    cg::grid_group grid,
    const bool active,
    const int packed,
    const int n_packed_threads,
    const int expansion_order,

    const CudaLincsBond* __restrict__ bonds,
    const int* __restrict__ neighbor_counts,
    const int* __restrict__ neighbor_indices,
    const double* __restrict__ matrix_a,
    const coord_t* __restrict__ directions,

    double* rhs_ping,
    double* rhs_pong,
    double* solution,

    const double* __restrict__ winv,
    coord_t* coords,

    int* failed) {
    double* rhs_in = rhs_ping;
    double* rhs_out = rhs_pong;

    for (int order = 0; order < expansion_order; order++) {
        if (active) {
            const int neighbor_count = neighbor_counts[packed];

            double next_value = 0;

            for (int slot = 0; slot < neighbor_count; slot++) {
                const int storage_idx = slot * n_packed_threads + packed;

                const int neighbor_idx = neighbor_indices[storage_idx];

                next_value += matrix_a[storage_idx] * rhs_in[neighbor_idx];
            }

            if (isfinite(next_value)) {
                rhs_out[packed] = next_value;
                solution[packed] += next_value;
            } else {
                rhs_out[packed] = 0;
                atomicCAS(failed, -1, packed);
            }
        }

        grid.sync();

        double* temporary = rhs_in;
        rhs_in = rhs_out;
        rhs_out = temporary;
    }

    if (active) {
        const auto& bond = bonds[packed];
        const coord_t direction = directions[packed];
        const double lambda = bond.mass_scale * solution[packed];

        if (isfinite(lambda)) {
            atomicAdd(&coords[bond.ai].x, -winv[bond.ai] * lambda * direction.x);
            atomicAdd(&coords[bond.ai].y, -winv[bond.ai] * lambda * direction.y);
            atomicAdd(&coords[bond.ai].z, -winv[bond.ai] * lambda * direction.z);

            atomicAdd(&coords[bond.aj].x, winv[bond.aj] * lambda * direction.x);
            atomicAdd(&coords[bond.aj].y, winv[bond.aj] * lambda * direction.y);
            atomicAdd(&coords[bond.aj].z, winv[bond.aj] * lambda * direction.z);
        } else {
            atomicCAS(failed, -1, packed);
        }
    }
    grid.sync();
}

__global__ void big_lincs_persistent_kernel(
    const int n_packed_threads,
    const int n_real_constraints,
    const int expansion_order,
    const int minimum_rotation_iterations,
    const int maximum_rotation_iterations,
    const double accuracy_tolerance,

    const CudaLincsBond* __restrict__ bonds,
    const int* __restrict__ neighbor_counts,
    const int* __restrict__ neighbor_indices,
    const double* __restrict__ matrix_a,
    const coord_t* __restrict__ directions,

    double* rhs_ping,
    double* rhs_pong,
    double* solution,

    const double* __restrict__ winv,
    coord_t* coords,

    int* failed,
    int* not_converged) {
    cg::grid_group grid = cg::this_grid();

    const int packed =
        blockIdx.x * blockDim.x + threadIdx.x;

    const bool active =
        packed < n_real_constraints;

    /*
     * Initial RHS and solution were produced by
     * big_lincs_prepare_direction_rhs_kernel.
     */
    big_lincs_solve_grid(
        grid,
        active,
        packed,
        n_packed_threads,
        expansion_order,
        bonds,
        neighbor_counts,
        neighbor_indices,
        matrix_a,
        directions,
        rhs_ping,
        rhs_pong,
        solution,
        winv,
        coords,
        failed);

    bool converged = false;

    for (int iteration = 0; iteration < maximum_rotation_iterations; iteration++) {
        /*
         * Construct rotational RHS and reset solution.
         */
        if (active) {
            const CudaLincsBond& bond = bonds[packed];

            const coord_t current_bond = coords[bond.ai] - coords[bond.aj];

            const double current_length2 = norm2(current_bond);

            const double target_length2 = bond.target * bond.target;

            const double radicand = 2.0 * target_length2 - current_length2;

            double rotation_rhs = 0;

            if (isfinite(radicand) && radicand > 0) {
                rotation_rhs = bond.mass_scale * (bond.target - sqrt(radicand));

                if (!isfinite(rotation_rhs)) {
                    rotation_rhs = 0;
                    atomicCAS(failed, -1, packed);
                }
            } else {
                atomicCAS(failed, -1, packed);
            }

            rhs_ping[packed] = rotation_rhs;
            solution[packed] = rotation_rhs;
        }

        /*
         * Every RHS entry must be initialized before expansion.
         */
        grid.sync();

        big_lincs_solve_grid(
            grid,
            active,
            packed,
            n_packed_threads,
            expansion_order,
            bonds,
            neighbor_counts,
            neighbor_indices,
            matrix_a,
            directions,
            rhs_ping,
            rhs_pong,
            solution,
            winv,
            coords,
            failed);

        if (iteration + 1 >= minimum_rotation_iterations) {
            /*
             * Reset the grid-wide reduction flag.
             */
            if (grid.thread_rank() == 0) {
                not_converged[0] = 0;
            }

            grid.sync();

            if (active) {
                const CudaLincsBond& bond = bonds[packed];

                const coord_t corrected_bond = coords[bond.ai] - coords[bond.aj];

                const double corrected_length = norm(corrected_bond);

                const double relative_error = fabs(corrected_length / bond.target - 1.0);

                if (!isfinite(corrected_length) || !isfinite(relative_error) || relative_error > accuracy_tolerance) {
                    atomicExch(not_converged, 1);
                }
            }

            grid.sync();

            /*
             * All threads observe the same flag and therefore
             * leave the loop uniformly.
             */
            if (not_converged[0] == 0) {
                converged = true;
                break;
            }
        }
    }

    /*
     * Report constraints that remain invalid after the maximum.
     */
    if (!converged && active) {
        const CudaLincsBond& bond = bonds[packed];

        const coord_t corrected_bond = coords[bond.ai] - coords[bond.aj];

        const double corrected_length = norm(corrected_bond);

        const double relative_error = fabs(corrected_length / bond.target - 1.0);

        if (!isfinite(corrected_length) || !isfinite(relative_error) || relative_error > accuracy_tolerance) {
            printf(">>> Lincs failed, i = %d, j = %d, d = %f, d0 = %f\n", bond.ai, bond.aj, corrected_length, bond.target);
            atomicCAS(failed, -1, packed);
        }
    }
}

__global__ void big_lincs_cooperative_solve_kernel(
    const int n_packed_threads,
    const int n_real_constraints,
    const int expansion_order,
    const int* neighbor_counts,
    const int* neighbor_indices,

    const double* matrix_a,
    double* rhs_in,
    double* rhs_out,
    double* solution,

    int* failed) {
    cg::grid_group grid = cg::this_grid();

    const int packed = blockIdx.x * blockDim.x + threadIdx.x;
    const bool active = packed < n_real_constraints;

    for (int i = 0; i < expansion_order; i++) {
        if (active) {
            const int neighbor_count = neighbor_counts[packed];
            double next_value = 0;

            for (int j = 0; j < neighbor_count; j++) {
                const int storage_idx = j * n_packed_threads + packed;
                const int neighbor_idx = neighbor_indices[storage_idx];
                next_value += matrix_a[storage_idx] * rhs_in[neighbor_idx];
            }

            if (isfinite(next_value)) {
                rhs_out[packed] = next_value;
                solution[packed] += next_value;
            } else {
                rhs_out[packed] = 0;
                atomicCAS(failed, -1, packed);
            }
        }

        grid.sync();
        double* temporary = rhs_in;
        rhs_in = rhs_out;
        rhs_out = temporary;
    }
}

__global__ void big_lincs_solve_kernel(
    const int n_packed_threads,
    const int n_real_constraints,
    const int* neighbor_counts,
    const int* neighbor_indices,

    double* matrix_a,
    double* rhs_in,
    double* rhs_out,
    double* solution,

    int* failed) {
    const int packed = blockIdx.x * blockDim.x + threadIdx.x;
    if (packed >= n_real_constraints) return;

    int neighbor_count = neighbor_counts[packed];
    double next_value = 0;
    for (int i = 0; i < neighbor_count; i++) {
        const int storage_idx = i * n_packed_threads + packed;
        const int neighbor_idx = neighbor_indices[storage_idx];

        next_value += matrix_a[storage_idx] * rhs_in[neighbor_idx];
    }
    rhs_out[packed] = next_value;
    solution[packed] += next_value;
}

__global__ void big_lincs_apply_correction_kernel(
    const int n_real_constraints,
    const CudaLincsBond* bonds,

    const coord_t* directions,
    const double* solution,

    const double* winv,
    coord_t* coords,
    int* failed) {
    const int packed = blockIdx.x * blockDim.x + threadIdx.x;
    if (packed >= n_real_constraints) return;

    const auto& bond = bonds[packed];
    const auto& direction = directions[packed];
    const double lambda = bond.mass_scale * solution[packed];
    if (isfinite(lambda)) {
        atomicAdd(&coords[bond.ai].x, -winv[bond.ai] * lambda * direction.x);
        atomicAdd(&coords[bond.ai].y, -winv[bond.ai] * lambda * direction.y);
        atomicAdd(&coords[bond.ai].z, -winv[bond.ai] * lambda * direction.z);
        atomicAdd(&coords[bond.aj].x, winv[bond.aj] * lambda * direction.x);
        atomicAdd(&coords[bond.aj].y, winv[bond.aj] * lambda * direction.y);
        atomicAdd(&coords[bond.aj].z, winv[bond.aj] * lambda * direction.z);
    } else {
        atomicCAS(failed, -1, packed);
    }
}

__global__ void big_lincs_rotation_rhs_kernel(
    const int n_real_constraints,
    const CudaLincsBond* bonds,

    double* rhs_ping,
    double* solution,

    coord_t* coords,
    int* failed) {
    const int packed = blockIdx.x * blockDim.x + threadIdx.x;
    if (packed >= n_real_constraints) return;

    const auto& bond = bonds[packed];
    const coord_t current_bond = coords[bond.ai] - coords[bond.aj];
    const double current_length2 = norm2(current_bond);
    const double target_length2 = bond.target * bond.target;
    const double radicand = 2.0 * target_length2 - current_length2;

    double rotation_rhs = 0;
    if (isfinite(radicand) && radicand > 0) {
        rotation_rhs = bond.mass_scale * (bond.target - sqrt(radicand));
        if (!isfinite(rotation_rhs)) {
            rotation_rhs = 0;
            atomicCAS(failed, -1, packed);
        }
    } else {
        rotation_rhs = 0;
        atomicCAS(failed, -1, packed);
    }
    rhs_ping[packed] = rotation_rhs;
    solution[packed] = rotation_rhs;
}

__global__ void big_lincs_validate_kernel(
    const int n_real_constraints,
    const double accuracy_tolerance,
    const CudaLincsBond* bonds,
    const coord_t* coords,
    int* failed) {
    const int packed = blockIdx.x * blockDim.x + threadIdx.x;
    if (packed >= n_real_constraints) return;

    const auto& bond = bonds[packed];
    const coord_t corrected_bond = coords[bond.ai] - coords[bond.aj];
    const double corrected_length = norm(corrected_bond);
    const double relative_error = fabs(corrected_length / bond.target - 1.0);

    if (!isfinite(corrected_length) || !isfinite(relative_error) || relative_error > accuracy_tolerance) {
        printf(">>> Lincs failed, i = %d, j = %d, d = %f, d0 = %f\n", bond.ai, bond.aj, corrected_length, bond.target);
        atomicCAS(failed, -1, packed);
    }
}

__global__ void big_lincs_check_convergence_kernel(
    const int n_real_constraints,
    const double accuracy_tolerance,
    const CudaLincsBond* bonds,
    const coord_t* coords,
    int* not_converged) {
    const int packed = blockIdx.x * blockDim.x + threadIdx.x;
    if (packed >= n_real_constraints) return;

    const auto& bond = bonds[packed];
    const coord_t corrected_bond = coords[bond.ai] - coords[bond.aj];
    const double corrected_length = norm(corrected_bond);
    const double relative_error = fabs(corrected_length / bond.target - 1.0);

    if (!isfinite(corrected_length) || !isfinite(relative_error) || relative_error > accuracy_tolerance) {
        atomicExch(not_converged, 1);
    }
}

}  // namespace

void CudaLincs::apply(Context& ctx, HostDeviceBuffer<coord_t>& xcoords) {
    coord_t* candidate = ctx.coords->gpu_data_p;
    const coord_t* reference = xcoords.gpu_data_p;
    apply_to(ctx, candidate, reference);
}

void CudaLincs::init_backend_for_big_group(
    Context& ctx,
    const std::vector<std::vector<int>>& bond_graph,
    const std::vector<int>& big_group_idx,
    const std::vector<std::vector<int>>& groups) {
    if (big_group_idx.empty()) return;

    const int bond_num = bond_graph.size();
    const auto* constraint_bonds = data_.constraint_bonds->cpu_data_p;
    const auto* winv = ctx.winv->cpu_data_p;

    /*
     * These are bond/constraint indices, not atom indices. Consecutive groups
     * of 256 constraints are owned by one CUDA block. A row's neighbors may
     * belong to any block and will be read from the global RHS buffers.
     */
    std::vector<int> bfs_order;
    std::vector<unsigned char> visited(bond_num, 0);

    for (int group_idx : big_group_idx) {
        const auto& group = groups[group_idx];
        if (group.empty()) continue;

        std::queue<int> pending;
        pending.push(group.front());
        visited[group.front()] = 1;

        while (!pending.empty()) {
            const int bond_idx = pending.front();
            pending.pop();
            bfs_order.push_back(bond_idx);

            for (int neighbor_idx : bond_graph[bond_idx]) {
                if (!visited[neighbor_idx]) {
                    visited[neighbor_idx] = 1;
                    pending.push(neighbor_idx);
                }
            }
        }
    }

    const int n_real_constraints = bfs_order.size();
    const int n_blocks = (n_real_constraints + BLOCK_THREAD - 1) / BLOCK_THREAD;
    const int n_packed_threads = n_blocks * BLOCK_THREAD;

    const CudaLincsBond dummy_bond = {-1, -1, 0, 0};
    std::vector<CudaLincsBond> bonds(n_packed_threads, dummy_bond);
    std::vector<int> neighbor_counts(n_packed_threads, 0);
    std::vector<int> original_to_packed(bond_num, -1);

    for (int packed = 0; packed < n_real_constraints; packed++) {
        const int original = bfs_order[packed];
        const ConstraintBond& bond = constraint_bonds[original];
        const int ai = bond.ai - 1;
        const int aj = bond.aj - 1;

        original_to_packed[original] = packed;
        bonds[packed] = {
            ai,
            aj,
            std::sqrt(bond.dist2),
            1.0 / std::sqrt(winv[ai] + winv[aj])};
        neighbor_counts[packed] = bond_graph[original].size();
    }

    int max_neighbors = 0;
    for (int packed = 0; packed < n_real_constraints; packed++) {
        max_neighbors = std::max(max_neighbors, neighbor_counts[packed]);
    }

    // Transposed ELLPACK keeps each neighbor slot coalesced across a warp.
    std::vector<int> neighbor_indices(max_neighbors * n_packed_threads, -1);
    std::vector<double> mass_factors(max_neighbors * n_packed_threads, 0);

    for (int packed = 0; packed < n_real_constraints; packed++) {
        const int original = bfs_order[packed];
        const ConstraintBond& bond = constraint_bonds[original];

        for (int slot = 0; slot < neighbor_counts[packed]; slot++) {
            const int neighbor_original = bond_graph[original][slot];
            const int neighbor_packed = original_to_packed[neighbor_original];
            const int storage_idx = slot * n_packed_threads + packed;

            if (neighbor_packed < 0) {
                throw std::runtime_error("CUDA LINCS large-group neighbor was not packed");
            }

            neighbor_indices[storage_idx] = neighbor_packed;
            const int shared_atom = get_same_atom_in_two_bond(bond, constraint_bonds[neighbor_original]);
            mass_factors[storage_idx] =
                -get_sign(shared_atom, bond) *
                get_sign(shared_atom, constraint_bonds[neighbor_original]) *
                winv[shared_atom - 1] *
                bonds[packed].mass_scale *
                bonds[neighbor_packed].mass_scale;
        }
    }

    int device = 0;
    check_cuda(cudaGetDevice(&device));

    int multiprocessor_count = 0;
    check_cuda(cudaDeviceGetAttribute(&multiprocessor_count, cudaDevAttrMultiProcessorCount, device));

    int active_blocks_per_multiprocessor = 0;
    check_cuda(cudaOccupancyMaxActiveBlocksPerMultiprocessor(&active_blocks_per_multiprocessor, big_lincs_cooperative_solve_kernel, BLOCK_THREAD, 0));

    const int cooperative_block_capacity = active_blocks_per_multiprocessor * multiprocessor_count;

    int cooperative_launch_supported = 0;
    check_cuda(cudaDeviceGetAttribute(&cooperative_launch_supported, cudaDevAttrCooperativeLaunch, device));
    const bool can_launch_cooperative_solve = n_blocks <= cooperative_block_capacity && cooperative_launch_supported != 0;

    int persistent_blocks_per_multiprocessor = 0;

    check_cuda(
        cudaOccupancyMaxActiveBlocksPerMultiprocessor(
            &persistent_blocks_per_multiprocessor,
            big_lincs_persistent_kernel,
            BLOCK_THREAD,
            0));

    const int persistent_block_capacity =
        persistent_blocks_per_multiprocessor *
        multiprocessor_count;

    CudaLincsBigData& data = lincs_big_data_;
    data.n_blocks = n_blocks;
    data.n_packed_threads = n_packed_threads;
    data.n_real_constraints = n_real_constraints;
    data.max_neighbors = max_neighbors;
    data.bonds = HostDeviceBuffer<CudaLincsBond>::from_vector(bonds, true);
    data.neighbor_counts = HostDeviceBuffer<int>::from_vector(neighbor_counts, true);
    data.neighbor_indices = HostDeviceBuffer<int>::from_vector(neighbor_indices, true);
    data.mass_factors = HostDeviceBuffer<double>::from_vector(mass_factors, true);
    data.directions = std::make_unique<HostDeviceBuffer<coord_t>>(n_packed_threads, false, true);
    data.matrix_a = std::make_unique<HostDeviceBuffer<double>>(mass_factors.size(), false, true);
    data.rhs_ping = std::make_unique<HostDeviceBuffer<double>>(n_packed_threads, false, true);
    data.rhs_pong = std::make_unique<HostDeviceBuffer<double>>(n_packed_threads, false, true);
    data.solution = std::make_unique<HostDeviceBuffer<double>>(n_packed_threads, false, true);
    data.failed_constraint = std::make_unique<HostDeviceBuffer<int>>(1, true, true);
    data.not_converged = std::make_unique<HostDeviceBuffer<int>>(1, true, true);
    data.can_launch_cooperative_solve = can_launch_cooperative_solve;
    data.can_launch_persistent = cooperative_launch_supported != 0 && n_blocks <= persistent_block_capacity;
}

void CudaLincs::init_backend(Context& ctx) {
    int bond_num = data_.n_constraints;
    auto* winv = ctx.winv->cpu_data_p;

    // 1. Get the correct bond group
    UnionFind uf(bond_num);

    auto* constraint_bonds = data_.constraint_bonds->cpu_data_p;

    std::vector<std::vector<int>> bond_graph(bond_num);

    for (int i = 0; i < bond_num; i++) {
        auto& bond_i = constraint_bonds[i];
        int a = bond_i.ai, b = bond_i.aj;
        for (int j = i + 1; j < bond_num; j++) {
            auto& bond_j = constraint_bonds[j];
            int a2 = bond_j.ai, b2 = bond_j.aj;
            if (a == a2 || a == b2 || b == a2 || b == b2) {
                uf.unite(i, j);
                bond_graph[i].push_back(j);
                bond_graph[j].push_back(i);
            }
        }
    }
    auto groups = uf.groups();
    // Build bonds
    std::vector<CudaLincsBond> bonds;
    std::vector<int> neighbor_counts;
    int group_num = groups.size();
    std::vector<int> f(group_num);
    std::iota(f.begin(), f.end(), 0);
    std::sort(f.begin(), f.end(), [&](int x, int y) {
        return groups[x].size() > groups[y].size();
    });

    int cur_bond_num = 0;
    std::vector<int> big_group_idx;

    CudaLincsBond dummy_bond = {-1, -1, 0, 0};

    std::vector<int> bond_thread_idx_map_to_idx;
    std::vector<int> bond_idx_map_to_thread_idx(bond_num, -1);
    for (int i = 0; i < group_num; i++) {
        int idx = f[i];
        int group_size = groups[idx].size();
        if (group_size > BLOCK_THREAD) {
            big_group_idx.push_back(idx);
            continue;
        }

        if (cur_bond_num + group_size > BLOCK_THREAD) {
            // add dummy first
            int left_size = BLOCK_THREAD - cur_bond_num;
            for (int j = 0; j < left_size; j++) {
                bond_thread_idx_map_to_idx.push_back(-1);
                bonds.push_back(dummy_bond);
                neighbor_counts.push_back(0);
            }
            cur_bond_num = 0;
        }

        // add in it
        for (int j = 0; j < group_size; j++) {
            int bond_idx = groups[idx][j];
            int ai = constraint_bonds[bond_idx].ai - 1, aj = constraint_bonds[bond_idx].aj - 1;
            bond_thread_idx_map_to_idx.push_back(bond_idx);
            bond_idx_map_to_thread_idx[bond_idx] = bonds.size();
            bonds.push_back({ai, aj, std::sqrt(constraint_bonds[bond_idx].dist2), 1.0 / sqrt(winv[ai] + winv[aj])});
            neighbor_counts.push_back(bond_graph[bond_idx].size());
        }
        cur_bond_num += group_size;
    }

    init_backend_for_big_group(ctx, bond_graph, big_group_idx, groups);

    // Pad the final block
    if (!bonds.empty() && cur_bond_num != 0) {
        for (int i = 0; i < BLOCK_THREAD - cur_bond_num; i++) {
            bond_thread_idx_map_to_idx.push_back(-1);
            bonds.push_back(dummy_bond);
            neighbor_counts.push_back(0);
        }
    }

    // Build neighbor_local_indices and mass factors
    int max_neighbors = 0;
    for (int count : neighbor_counts) {
        max_neighbors = std::max(max_neighbors, count);
    }

    int total_bond_num = bonds.size();
    std::vector<int> neighbor_local_indices(max_neighbors * total_bond_num, -1);
    std::vector<double> mass_factors(max_neighbors * total_bond_num, 0);

    for (int i = 0; i < total_bond_num; i++) {
        int bond_idx = bond_thread_idx_map_to_idx[i];
        if (bond_idx == -1) continue;  // dummy

        for (int j = 0; j < bond_graph[bond_idx].size(); j++) {
            int bond_idx2 = bond_graph[bond_idx][j];
            int bond_idx2_thread_idx = bond_idx_map_to_thread_idx[bond_idx2];

            int storage_idx = j * total_bond_num + i;
            neighbor_local_indices[storage_idx] = bond_idx2_thread_idx % BLOCK_THREAD;

            int atom_idx = get_same_atom_in_two_bond(constraint_bonds[bond_idx], constraint_bonds[bond_idx2]);
            mass_factors[storage_idx] = -get_sign(atom_idx, constraint_bonds[bond_idx]) * get_sign(atom_idx, constraint_bonds[bond_idx2]) * winv[atom_idx - 1] * bonds[i].mass_scale * bonds[bond_idx2_thread_idx].mass_scale;
        }
    }

    int n_real_constraints = 0;
    for (int original : bond_thread_idx_map_to_idx) {
        if (original >= 0) n_real_constraints++;
    }

    int n_packed_threads = bonds.size();
    lincs_small_data_.n_blocks = n_packed_threads / BLOCK_THREAD;  // Must n_packed_threads % BLOCK_THREAD = 0
    lincs_small_data_.n_packed_threads = n_packed_threads;
    lincs_small_data_.n_real_constraints = n_real_constraints;
    lincs_small_data_.max_neighbors = max_neighbors;
    lincs_small_data_.bonds = HostDeviceBuffer<CudaLincsBond>::from_vector(bonds, true);
    lincs_small_data_.neighbor_counts = HostDeviceBuffer<int>::from_vector(neighbor_counts, true);
    lincs_small_data_.neighbor_local_indices = HostDeviceBuffer<int>::from_vector(neighbor_local_indices, true);
    lincs_small_data_.mass_factors = HostDeviceBuffer<double>::from_vector(mass_factors, true);
    lincs_small_data_.matrix_a = std::make_unique<HostDeviceBuffer<double>>(mass_factors.size(), false, true);
    lincs_small_data_.failed_constraint = std::make_unique<HostDeviceBuffer<int>>(1, true, true);
}

void CudaLincs::initial_constraint(Context& ctx) {
    // Same idea as the CudaShake
    HostDeviceBuffer<coord_t> xcoords(ctx.n_atoms, true, true);
    auto* d_coords = ctx.coords->gpu_data_p;
    auto* d_xcoords = xcoords.gpu_data_p;

    check_cuda(cudaMemcpy(d_xcoords, d_coords, ctx.n_atoms * sizeof(coord_t), cudaMemcpyDeviceToDevice));
    apply_to(ctx, d_coords, d_xcoords);

    ctx.coords->download();

    auto* coords = ctx.coords->cpu_data_p;
    auto* velocities = ctx.velocities->cpu_data_p;
    auto* host_xcoords = xcoords.cpu_data_p;
    for (int i = 0; i < ctx.n_atoms; i++) {
        const double dt = ctx.dt;
        host_xcoords[i].x = coords[i].x - dt * velocities[i].x;
        host_xcoords[i].y = coords[i].y - dt * velocities[i].y;
        host_xcoords[i].z = coords[i].z - dt * velocities[i].z;
    }

    xcoords.upload();

    apply_to(ctx, d_xcoords, d_coords);

    xcoords.download();

    for (int i = 0; i < ctx.n_atoms; i++) {
        const double dt = ctx.dt;
        velocities[i].x = (coords[i].x - host_xcoords[i].x) / dt;
        velocities[i].y = (coords[i].y - host_xcoords[i].y) / dt;
        velocities[i].z = (coords[i].z - host_xcoords[i].z) / dt;
    }
    ctx.velocities->upload();
}

void CudaLincs::apply_to(Context& ctx, coord_t* coords, const coord_t* xcoords) {
    CudaLincsSmallData& data = lincs_small_data_;

    if (data.n_blocks > 0) {
        data.failed_constraint->cpu_data_p[0] = -1;
        data.failed_constraint->upload();
        small_lincs_kernel<<<data.n_blocks, BLOCK_THREAD>>>(
            data.n_packed_threads,
            lincs_settings_.expansion_order,
            lincs_settings_.minimum_rotation_iterations,
            lincs_settings_.maximum_rotation_iterations,
            lincs_settings_.accuracy_tolerance,
            data.bonds->gpu_data_p,
            data.neighbor_counts->gpu_data_p,
            data.neighbor_local_indices->gpu_data_p,
            data.mass_factors->gpu_data_p,
            data.matrix_a->gpu_data_p,
            ctx.winv->gpu_data_p,
            coords,
            xcoords,
            data.failed_constraint->gpu_data_p);
        check_cuda(cudaGetLastError());
        data.failed_constraint->download();
        const int failed = data.failed_constraint->cpu_data_p[0];
        if (failed != -1) {
            std::fflush(stdout);
            std::exit(EXIT_FAILURE);
        }
    }

    auto& big_lincs_data = lincs_big_data_;

    if (big_lincs_data.n_blocks > 0) {
        big_lincs_data.failed_constraint->cpu_data_p[0] = -1;
        big_lincs_data.failed_constraint->upload();

        big_lincs_prepare_direction_rhs_kernel<<<big_lincs_data.n_blocks, BLOCK_THREAD>>>(
            big_lincs_data.n_real_constraints,
            big_lincs_data.bonds->gpu_data_p,
            big_lincs_data.directions->gpu_data_p,
            big_lincs_data.rhs_ping->gpu_data_p,
            big_lincs_data.solution->gpu_data_p,
            coords,
            xcoords,
            big_lincs_data.failed_constraint->gpu_data_p);

        check_cuda(cudaGetLastError());

        big_lincs_prepare_matrix_kernel<<<big_lincs_data.n_blocks, BLOCK_THREAD>>>(
            big_lincs_data.n_packed_threads,
            big_lincs_data.n_real_constraints,
            big_lincs_data.neighbor_counts->gpu_data_p,
            big_lincs_data.neighbor_indices->gpu_data_p,
            big_lincs_data.mass_factors->gpu_data_p,
            big_lincs_data.directions->gpu_data_p,
            big_lincs_data.matrix_a->gpu_data_p,
            big_lincs_data.failed_constraint->gpu_data_p);

        check_cuda(cudaGetLastError());

        if (big_lincs_data.can_launch_persistent) {
            int n_packed_threads =
                big_lincs_data.n_packed_threads;

            int n_real_constraints =
                big_lincs_data.n_real_constraints;

            int expansion_order =
                lincs_settings_.expansion_order;

            int minimum_rotation_iterations =
                lincs_settings_.minimum_rotation_iterations;

            int maximum_rotation_iterations =
                lincs_settings_.maximum_rotation_iterations;

            double accuracy_tolerance =
                lincs_settings_.accuracy_tolerance;

            const CudaLincsBond* bonds =
                big_lincs_data.bonds->gpu_data_p;

            const int* neighbor_counts =
                big_lincs_data.neighbor_counts->gpu_data_p;

            const int* neighbor_indices =
                big_lincs_data.neighbor_indices->gpu_data_p;

            const double* matrix_a =
                big_lincs_data.matrix_a->gpu_data_p;

            const coord_t* directions =
                big_lincs_data.directions->gpu_data_p;

            double* rhs_ping =
                big_lincs_data.rhs_ping->gpu_data_p;

            double* rhs_pong =
                big_lincs_data.rhs_pong->gpu_data_p;

            double* solution =
                big_lincs_data.solution->gpu_data_p;

            const double* winv =
                ctx.winv->gpu_data_p;

            int* failed =
                big_lincs_data.failed_constraint->gpu_data_p;

            int* not_converged =
                big_lincs_data.not_converged->gpu_data_p;

            void* kernel_arguments[] = {
                &n_packed_threads,
                &n_real_constraints,
                &expansion_order,
                &minimum_rotation_iterations,
                &maximum_rotation_iterations,
                &accuracy_tolerance,
                &bonds,
                &neighbor_counts,
                &neighbor_indices,
                &matrix_a,
                &directions,
                &rhs_ping,
                &rhs_pong,
                &solution,
                &winv,
                &coords,
                &failed,
                &not_converged,
            };

            check_cuda(cudaLaunchCooperativeKernel(
                reinterpret_cast<void*>(
                    big_lincs_persistent_kernel),
                dim3(big_lincs_data.n_blocks),
                dim3(BLOCK_THREAD),
                kernel_arguments,
                0,
                nullptr));

            check_cuda(cudaGetLastError());

            big_lincs_data.failed_constraint->download();

            const int failed_index =
                big_lincs_data
                    .failed_constraint
                    ->cpu_data_p[0];

            if (failed_index != -1) {
                std::fflush(stdout);
                std::exit(EXIT_FAILURE);
            }
        } else {
            auto solve_expansion = [&]() {
                if (big_lincs_data.can_launch_cooperative_solve) {
                    int n_packed_threads = big_lincs_data.n_packed_threads;

                    int n_real_constraints = big_lincs_data.n_real_constraints;

                    int expansion_order = lincs_settings_.expansion_order;

                    const int* neighbor_counts = big_lincs_data.neighbor_counts->gpu_data_p;

                    const int* neighbor_indices = big_lincs_data.neighbor_indices->gpu_data_p;

                    const double* matrix_a = big_lincs_data.matrix_a->gpu_data_p;

                    double* rhs_ping = big_lincs_data.rhs_ping->gpu_data_p;

                    double* rhs_pong = big_lincs_data.rhs_pong->gpu_data_p;

                    double* solution = big_lincs_data.solution->gpu_data_p;

                    int* failed = big_lincs_data.failed_constraint->gpu_data_p;

                    void* kernel_arguments[] = {
                        &n_packed_threads,
                        &n_real_constraints,
                        &expansion_order,
                        &neighbor_counts,
                        &neighbor_indices,
                        &matrix_a,
                        &rhs_ping,
                        &rhs_pong,
                        &solution,
                        &failed,
                    };

                    check_cuda(cudaLaunchCooperativeKernel(
                        reinterpret_cast<void*>(big_lincs_cooperative_solve_kernel),
                        dim3(big_lincs_data.n_blocks),
                        dim3(BLOCK_THREAD),
                        kernel_arguments,
                        0,
                        nullptr));

                    check_cuda(cudaGetLastError());

                } else {
                    double* rhs_in = big_lincs_data.rhs_ping->gpu_data_p;
                    double* rhs_out = big_lincs_data.rhs_pong->gpu_data_p;

                    for (int order = 0; order < lincs_settings_.expansion_order; order++) {
                        big_lincs_solve_kernel<<<big_lincs_data.n_blocks, BLOCK_THREAD>>>(
                            big_lincs_data.n_packed_threads,
                            big_lincs_data.n_real_constraints,
                            big_lincs_data.neighbor_counts->gpu_data_p,
                            big_lincs_data.neighbor_indices->gpu_data_p,
                            big_lincs_data.matrix_a->gpu_data_p,
                            rhs_in,
                            rhs_out,
                            big_lincs_data.solution->gpu_data_p,
                            big_lincs_data.failed_constraint->gpu_data_p);
                        check_cuda(cudaGetLastError());
                        std::swap(rhs_in, rhs_out);
                    }
                }
            };

            auto apply_correction = [&]() {
                big_lincs_apply_correction_kernel<<<big_lincs_data.n_blocks, BLOCK_THREAD>>>(
                    big_lincs_data.n_real_constraints,
                    big_lincs_data.bonds->gpu_data_p,
                    big_lincs_data.directions->gpu_data_p,
                    big_lincs_data.solution->gpu_data_p,
                    ctx.winv->gpu_data_p,
                    coords,
                    big_lincs_data.failed_constraint->gpu_data_p);
                check_cuda(cudaGetLastError());
            };

            solve_expansion();
            apply_correction();

            bool converged = false;
            for (int iteration = 0; iteration < lincs_settings_.maximum_rotation_iterations; iteration++) {
                big_lincs_rotation_rhs_kernel<<<big_lincs_data.n_blocks, BLOCK_THREAD>>>(
                    big_lincs_data.n_real_constraints,
                    big_lincs_data.bonds->gpu_data_p,
                    big_lincs_data.rhs_ping->gpu_data_p,
                    big_lincs_data.solution->gpu_data_p,
                    coords,
                    big_lincs_data.failed_constraint->gpu_data_p);
                check_cuda(cudaGetLastError());

                solve_expansion();
                apply_correction();

                if (iteration + 1 >= lincs_settings_.minimum_rotation_iterations) {
                    big_lincs_data.not_converged->cpu_data_p[0] = 0;
                    big_lincs_data.not_converged->upload();

                    big_lincs_check_convergence_kernel<<<big_lincs_data.n_blocks, BLOCK_THREAD>>>(
                        big_lincs_data.n_real_constraints,
                        lincs_settings_.accuracy_tolerance,
                        big_lincs_data.bonds->gpu_data_p,
                        coords,
                        big_lincs_data.not_converged->gpu_data_p);
                    check_cuda(cudaGetLastError());

                    big_lincs_data.not_converged->download();
                    if (big_lincs_data.not_converged->cpu_data_p[0] == 0) {
                        converged = true;
                        break;
                    }
                }
            }

            if (!converged) {
                big_lincs_validate_kernel<<<big_lincs_data.n_blocks, BLOCK_THREAD>>>(
                    big_lincs_data.n_real_constraints,
                    lincs_settings_.accuracy_tolerance,
                    big_lincs_data.bonds->gpu_data_p,
                    coords,
                    big_lincs_data.failed_constraint->gpu_data_p);
                check_cuda(cudaGetLastError());
            }

            big_lincs_data.failed_constraint->download();
            const int failed = big_lincs_data.failed_constraint->cpu_data_p[0];
            if (failed != -1) {
                std::fflush(stdout);
                std::exit(EXIT_FAILURE);
            }
        }
    }
}

void CudaLincs::cleanup() {
}
