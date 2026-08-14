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
    const int begin,
    const int end,
    const int n_constraints,  // total constraints
    const int expansion_order,
    const int* neighbor_counts,  // neighbor count
    const int* neighbor_indices,
    const CudaLincsBond* bonds,
    const double* matrix_a,
    const double* winv,
    const coord_t* directions,

    double* q0,
    double* q1,
    double* solution,
    coord_t* coords,
    int* failed) {
    /*
     * Neumann expansion:
     *
     * solution =
     *     rhs + A*rhs + A^2*rhs + ... + A^order*rhs
     */

    for (int constraint = begin + threadIdx.x; constraint < end; constraint += blockDim.x) {
        solution[constraint] = q0[constraint];
    }
    __syncthreads();

    double* input = q0;
    double* output = q1;

    for (int i = 0; i < expansion_order; i++) {
        for (int constraint = begin + threadIdx.x; constraint < end; constraint += blockDim.x) {
            double next = 0;

            const int count = neighbor_counts[constraint];
            for (int j = 0; j < count; j++) {
                const int storage_idx = j * n_constraints + constraint;

                const int neighbor = neighbor_indices[storage_idx];

                next += matrix_a[storage_idx] * input[neighbor];
            }
            output[constraint] = next;
            solution[constraint] += next;
        }
        __syncthreads();

        double* temp = input;
        input = output;
        output = temp;
    }

    /*
     * Apply the initial coordinate correction:
     *
     * lambda = mass_scale * solution
     *
     * ai -= winv[ai] * lambda * direction
     * aj += winv[aj] * lambda * direction
     */

    for (int constraint = begin + threadIdx.x; constraint < end; constraint += blockDim.x) {
        const CudaLincsBond& bond = bonds[constraint];
        const double lambda = bond.mass_scale * solution[constraint];
        const coord_t direction = directions[constraint];

        if (isfinite(lambda)) {
            atomicAdd(&coords[bond.ai].x, -winv[bond.ai] * lambda * direction.x);
            atomicAdd(&coords[bond.ai].y, -winv[bond.ai] * lambda * direction.y);
            atomicAdd(&coords[bond.ai].z, -winv[bond.ai] * lambda * direction.z);
            atomicAdd(&coords[bond.aj].x, winv[bond.aj] * lambda * direction.x);
            atomicAdd(&coords[bond.aj].y, winv[bond.aj] * lambda * direction.y);
            atomicAdd(&coords[bond.aj].z, winv[bond.aj] * lambda * direction.z);
        } else {
            atomicCAS(failed, -1, threadIdx.x);
        }
    }
}

__global__ void lincs_kernel(
    const int n_constraints,
    const int expansion_order,
    const int minimum_rotation_iterations,
    const int maximum_rotation_iterations,
    const double accuracy_tolerance,
    const int* work_unit_offsets,
    const CudaLincsBond* bonds,
    const int* neighbor_counts,
    const int* neighbor_indices,
    const double* mass_factors,

    coord_t* directions,
    double* matrix_a,
    double* q0,
    double* q1,
    double* solution,

    const double* winv,
    coord_t* coords,
    const coord_t* reference,
    int* failed

) {
    const int local = threadIdx.x;
    const int block_idx = blockIdx.x;
    const int begin = work_unit_offsets[block_idx];
    const int end = work_unit_offsets[block_idx + 1];
    __shared__ int block_ready;

    for (int constraint = begin + threadIdx.x; constraint < end; constraint += blockDim.x) {
        const CudaLincsBond& bond = bonds[constraint];
        const coord_t old_bond = reference[bond.ai] - reference[bond.aj];
        const double old_d2 = norm2(old_bond);

        coord_t direction = {0, 0, 0};
        if (isfinite(old_d2) && old_d2 > 0) {
            const double inverse_length = rsqrt(old_d2);
            direction = old_bond * inverse_length;
        } else {
            atomicCAS(failed, -1, constraint);
        }
        directions[constraint] = direction;
    }
    __syncthreads();

    // Construct matrix-A
    for (int constraint = begin + threadIdx.x; constraint < end; constraint += blockDim.x) {
        const CudaLincsBond& bond = bonds[constraint];
        const coord_t direction = directions[constraint];

        for (int j = 0; j < neighbor_counts[constraint]; j++) {
            const int storage_idx = j * n_constraints + constraint;
            const int neighbor = neighbor_indices[storage_idx];
            matrix_a[storage_idx] = mass_factors[storage_idx] * dot(direction, directions[neighbor]);
        }

        const coord_t current_bond = coords[bond.ai] - coords[bond.aj];
        double rhs = bond.mass_scale * (dot(direction, current_bond) - bond.target);
        if (!isfinite(rhs)) {
            rhs = 0;
            atomicCAS(failed, -1, constraint);
        }
        q0[constraint] = rhs;
    }
    solve(begin, end, n_constraints, expansion_order, neighbor_counts, neighbor_indices, bonds, matrix_a, winv, directions, q0, q1, solution, coords, failed);
    __syncthreads();
    /*
     * Rotational/centripetal corrections.
     */
    for (int iteration = 0; iteration < maximum_rotation_iterations; iteration++) {
        for (int constraint = begin + threadIdx.x; constraint < end; constraint += blockDim.x) {
            double rotation_rhs = 0;
            const CudaLincsBond& bond = bonds[constraint];
            const coord_t current_bond = coords[bond.ai] - coords[bond.aj];
            const double current_length2 = norm2(current_bond);

            const double target_length2 = bond.target * bond.target;

            const double radicand = 2.0 * target_length2 - current_length2;
            if (isfinite(radicand) && radicand > 0) {
                rotation_rhs = bond.mass_scale * (bond.target - sqrt(radicand));
                if (!isfinite(rotation_rhs)) {
                    rotation_rhs = 0;
                    atomicCAS(failed, -1, constraint);
                }
            } else {
                rotation_rhs = 0;
                atomicCAS(failed, -1, constraint);
            }
            q0[constraint] = rotation_rhs;
        }
        solve(begin, end, n_constraints, expansion_order, neighbor_counts, neighbor_indices, bonds, matrix_a, winv, directions, q0, q1, solution, coords, failed);

        if (local == 0) {
            block_ready = 1;
        }
        __syncthreads();

        for (int constraint = begin + threadIdx.x; constraint < end; constraint += blockDim.x) {
            const CudaLincsBond& bond = bonds[constraint];
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
    for (int constraint = begin + threadIdx.x; constraint < end; constraint += blockDim.x) {
        const CudaLincsBond& bond = bonds[constraint];
        const coord_t corrected_bond = coords[bond.ai] - coords[bond.aj];
        const double corrected_length = norm(corrected_bond);

        const double relative_error = fabs(corrected_length / bond.target - 1.0);
        if (!isfinite(corrected_length) || !isfinite(relative_error) || relative_error > accuracy_tolerance) {
            printf(">>> Lincs failed, i = %d, j = %d, d = %f, d0 = %f\n", bond.ai, bond.aj, corrected_length, bond.target);
            atomicCAS(failed, -1, constraint);
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
    int group_num = groups.size();
    std::vector<int> f(group_num);
    std::iota(f.begin(), f.end(), 0);
    std::sort(f.begin(), f.end(), [&](int x, int y) {
        return groups[x].size() > groups[y].size();
    });

    // Build bonds
    std::vector<CudaLincsBond> bonds;
    std::vector<int> bond_idx_to_thread_idx(bond_num);
    for (int i = 0; i < group_num; i++) {
        int idx = f[i];
        int group_size = groups[idx].size();
        for (int j = 0; j < group_size; j++) {
            int bond_idx = groups[idx][j];
            int ai = constraint_bonds[bond_idx].ai - 1, aj = constraint_bonds[bond_idx].aj - 1;
            bond_idx_to_thread_idx[bond_idx] = bonds.size();
            bonds.push_back({ai, aj, std::sqrt(constraint_bonds[bond_idx].dist2), 1.0 / sqrt(winv[ai] + winv[aj])});
        }
    }

    // Build work_unit_offsets
    int cur_bond_num = 0;
    std::vector<int> work_unit_offsets;
    work_unit_offsets.push_back(0);
    int max_work_unit = 0;

    auto close_work_unit = [&]() {
        if (cur_bond_num == 0) {
            return;
        }

        max_work_unit = std::max(max_work_unit, cur_bond_num);
        work_unit_offsets.push_back(work_unit_offsets.back() + cur_bond_num);
        cur_bond_num = 0;
    };

    for (int i = 0; i < group_num; i++) {
        int idx = f[i];
        int group_size = groups[idx].size();

        if (group_size > BLOCK_THREAD) {
            close_work_unit();
            cur_bond_num = group_size;
            close_work_unit();
        } else {
            if (cur_bond_num + group_size > BLOCK_THREAD) {
                close_work_unit();
            }
            cur_bond_num += group_size;
        }
    }
    close_work_unit();

    // Build neighbor_counts
    std::vector<int> neighbor_counts;
    for (int i = 0; i < group_num; i++) {
        int idx = f[i];
        int group_size = groups[idx].size();
        for (int j = 0; j < group_size; j++) {
            int bond_idx = groups[idx][j];
            neighbor_counts.push_back(bond_graph[bond_idx].size());
        }
    }

    // Build neighbor_local_indices and mass factors
    int max_neighbors = 0;
    for (int count : neighbor_counts) {
        max_neighbors = std::max(max_neighbors, count);
    }

    // Build neighbor_indices
    // Build mass factors
    std::vector<int> neighbor_indices(max_neighbors * bond_num, -1);
    std::vector<double> mass_factors(max_neighbors * bond_num, 0);
    int cur_thread_idx = 0;
    for (int i = 0; i < group_num; i++) {
        int idx = f[i];
        int group_size = groups[idx].size();
        for (int j = 0; j < group_size; j++) {
            int bond_idx = groups[idx][j];
            for (int k = 0; k < bond_graph[bond_idx].size(); k++) {
                int bond_idx2 = bond_graph[bond_idx][k];
                int bond_idx2_thread_idx = bond_idx_to_thread_idx[bond_idx2];
                int storage_idx = k * bond_num + cur_thread_idx;
                neighbor_indices[storage_idx] = bond_idx2_thread_idx;
                int atom_idx = get_same_atom_in_two_bond(constraint_bonds[bond_idx], constraint_bonds[bond_idx2]);
                mass_factors[storage_idx] = -get_sign(atom_idx, constraint_bonds[bond_idx]) * get_sign(atom_idx, constraint_bonds[bond_idx2]) * winv[atom_idx - 1] * bonds[cur_thread_idx].mass_scale * bonds[bond_idx2_thread_idx].mass_scale;
            }

            cur_thread_idx++;
        }
    }

    lincs_data_.n_blocks = work_unit_offsets.size() - 1;
    lincs_data_.n_constraints = bond_num;
    lincs_data_.max_neighbors = max_neighbors;
    lincs_data_.max_work_unit_size = max_work_unit;
    lincs_data_.work_unit_offsets = HostDeviceBuffer<int>::from_vector(work_unit_offsets, true);
    lincs_data_.bonds = HostDeviceBuffer<CudaLincsBond>::from_vector(bonds, true);
    lincs_data_.neighbor_counts = HostDeviceBuffer<int>::from_vector(neighbor_counts, true);
    lincs_data_.neighbor_indices = HostDeviceBuffer<int>::from_vector(neighbor_indices, true);
    lincs_data_.mass_factors = HostDeviceBuffer<double>::from_vector(mass_factors, true);

    lincs_data_.directions = std::make_unique<HostDeviceBuffer<coord_t>>(bond_num, false, true);
    lincs_data_.matrix_a = std::make_unique<HostDeviceBuffer<double>>(mass_factors.size(), false, true);
    lincs_data_.q0 = std::make_unique<HostDeviceBuffer<double>>(bond_num, false, true);
    lincs_data_.q1 = std::make_unique<HostDeviceBuffer<double>>(bond_num, false, true);
    lincs_data_.solution = std::make_unique<HostDeviceBuffer<double>>(bond_num, false, true);
    lincs_data_.failed_constraint = std::make_unique<HostDeviceBuffer<int>>(1, true, true);
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
    CudaLincsData& data = lincs_data_;

    if (data.n_blocks == 0) return;

    data.failed_constraint->cpu_data_p[0] = -1;
    data.failed_constraint->upload();
    lincs_kernel<<<data.n_blocks, BLOCK_THREAD>>>(
        data.n_constraints,
        lincs_settings_.expansion_order,
        lincs_settings_.minimum_rotation_iterations,
        lincs_settings_.maximum_rotation_iterations,
        lincs_settings_.accuracy_tolerance,
        data.work_unit_offsets->gpu_data_p,
        data.bonds->gpu_data_p,
        data.neighbor_counts->gpu_data_p,
        data.neighbor_indices->gpu_data_p,
        data.mass_factors->gpu_data_p,

        data.directions->gpu_data_p,
        data.matrix_a->gpu_data_p,
        data.q0->gpu_data_p,
        data.q1->gpu_data_p,
        data.solution->gpu_data_p,
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