#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <stdexcept>
#include <vector>

#include "constants.h"
#include "cuda_runtime_utility.h"
#include "cuda_serial_shake.cuh"

namespace {

constexpr int kSerialShakeThreads = 128;

__global__ void serial_shake_by_molecule_kernel(
    int n_molecules,
    int n_atoms,
    int n_constraints,
    const int* molecule_constraint_offsets,
    const ShakeBond* shake_bonds,
    coord_t* coords,
    const coord_t* xcoords,
    const double* winv,
    std::uint8_t* ready,
    int* failed) {
    const int mol = blockIdx.x * blockDim.x + threadIdx.x;
    if (mol >= n_molecules) return;
    const int replica = blockIdx.y;
    const int atom_offset = replica * n_atoms;
    const int ready_offset = replica * n_constraints;
    coords += atom_offset;
    xcoords += atom_offset;
    ready += ready_offset;
    failed += replica;

    const int first = molecule_constraint_offsets[mol];
    const int end = molecule_constraint_offsets[mol + 1];
    if (first == end) return;

    bool converged = false;
    int try_times = 0;
    do {
        converged = true;
        for (int shake = first; shake < end; shake++) {
            if (ready[shake]) continue;

            const ShakeBond shake_bond = shake_bonds[shake];
            const int ai = shake_bond.ai - 1;
            const int aj = shake_bond.aj - 1;

            const double xij_x = coords[ai].x - coords[aj].x;
            const double xij_y = coords[ai].y - coords[aj].y;
            const double xij_z = coords[ai].z - coords[aj].z;
            const double current_dist2 = xij_x * xij_x + xij_y * xij_y + xij_z * xij_z;
            const double diff = shake_bond.dist2 - current_dist2;

            if (fabs(diff) < shake_tol * shake_bond.dist2) {
                ready[shake] = 1;
            } else {
                converged = false;
            }

            // Match CpuShake/Q exactly: even a constraint that becomes ready
            // on this visit receives its final correction.
            const double xxij_x = xcoords[ai].x - xcoords[aj].x;
            const double xxij_y = xcoords[ai].y - xcoords[aj].y;
            const double xxij_z = xcoords[ai].z - xcoords[aj].z;
            const double scp = xij_x * xxij_x + xij_y * xxij_y + xij_z * xxij_z;
            const double winv_ai = winv[ai];
            const double winv_aj = winv[aj];
            const double corr = diff / (2.0 * scp * (winv_ai + winv_aj));

            coords[ai].x += xxij_x * corr * winv_ai;
            coords[ai].y += xxij_y * corr * winv_ai;
            coords[ai].z += xxij_z * corr * winv_ai;

            coords[aj].x -= xxij_x * corr * winv_aj;
            coords[aj].y -= xxij_y * corr * winv_aj;
            coords[aj].z -= xxij_z * corr * winv_aj;
        }
        try_times++;
    } while (try_times < shake_max_iter && !converged);

    if (!converged) {
        atomicExch(failed, 1);
        for (int shake = first; shake < end; shake++) {
            const int ai = shake_bonds[shake].ai - 1;
            const int aj = shake_bonds[shake].aj - 1;
            const double dx = xcoords[ai].x - xcoords[aj].x;
            const double dy = xcoords[ai].y - xcoords[aj].y;
            const double dz = xcoords[ai].z - xcoords[aj].z;
            printf(">>> Shake failed, replica = %d, i = %d, j = %d, d = %f, d0 = %f\n",
                   replica,
                   ai,
                   aj,
                   sqrt(dx * dx + dy * dy + dz * dz),
                   sqrt(shake_bonds[shake].dist2));
        }
    }
}

}  // namespace

void CudaSerialConstraintSolver::init(
    const std::vector<int>& molecule_constraint_counts,
    const std::vector<ShakeBond>& shake_bonds, int n_atoms, int max_replicates) {
    n_molecules_ = static_cast<int>(molecule_constraint_counts.size());
    n_constraints_ = static_cast<int>(shake_bonds.size());
    n_atoms_ = n_atoms;
    max_replicates_ = max_replicates;

    std::vector<int> offsets(n_molecules_ + 1, 0);
    for (int mol = 0; mol < n_molecules_; mol++) {
        offsets[mol + 1] = offsets[mol] + molecule_constraint_counts[mol];
    }
    if (offsets.back() != n_constraints_) {
        throw std::runtime_error("CUDA serial SHAKE molecule counts do not match the constraint list.");
    }

    molecule_constraint_offsets_ = HostDeviceBuffer<int>::from_vector(offsets, true);
    shake_bonds_ = HostDeviceBuffer<ShakeBond>::from_vector(shake_bonds, true);
    ready_ = std::make_unique<HostDeviceBuffer<std::uint8_t>>(n_constraints_ * max_replicates_, false, true);
    failed_ = std::make_unique<HostDeviceBuffer<int>>(max_replicates_, true, true);
}

void CudaSerialConstraintSolver::apply(
    coord_t* d_coords,
    const coord_t* d_xcoords,
    const double* d_winv, int n_replicates) {
    if (!enabled()) return;

    ready_->zero();
    failed_->zero();

    const int blocks = (n_molecules_ + kSerialShakeThreads - 1) / kSerialShakeThreads;
    serial_shake_by_molecule_kernel<<<dim3(blocks, n_replicates), kSerialShakeThreads>>>(
        n_molecules_, n_atoms_, n_constraints_,
        molecule_constraint_offsets_->gpu_data_p,
        shake_bonds_->gpu_data_p,
        d_coords,
        d_xcoords,
        d_winv,
        ready_->gpu_data_p,
        failed_->gpu_data_p);
    check_cuda(cudaGetLastError());

    // This synchronizes the SHAKE kernel, but transfers only one integer.  It
    // retains CpuShake's fail-fast behavior without moving coordinate arrays.
    failed_->download();
    for (int replica = 0; replica < n_replicates; replica++) {
        if (failed_->cpu_data_p[replica]) {
            std::exit(EXIT_FAILURE);
        }
    }
}
