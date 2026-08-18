#include "cuda_lincs.cuh"
#include "geometry.h"
#include "union_find.h"
#include "vector"

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

}  // namespace

void CudaLincs::apply(Context& ctx, HostDeviceBuffer<coord_t>& xcoords) {
    coord_t* candidate = ctx.coords->gpu_data_p;
    const coord_t* reference = xcoords.gpu_data_p;
    apply_to(ctx, candidate, reference);
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

    if (data.n_blocks == 0) return;

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

void CudaLincs::cleanup() {
}