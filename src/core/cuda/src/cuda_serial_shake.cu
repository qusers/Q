#include "cuda_serial_shake.cuh"

#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <stdexcept>
#include <vector>

#include "constants.h"
#include "cuda_runtime_utility.h"

namespace {

constexpr int kSerialShakeThreads = 128;

__global__ void serial_shake_by_molecule_kernel(
    int n_molecules,
    const int* molecule_constraint_offsets,
    const ShakeBond* shake_bonds,
    coord_t* coords,
    const coord_t* xcoords,
    const double* winv,
    std::uint8_t* ready,
    int* failed) {
    const int mol = blockIdx.x * blockDim.x + threadIdx.x;
    if (mol >= n_molecules) return;

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
            printf(">>> Shake failed, i = %d, j = %d, d = %f, d0 = %f\n",
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
    const std::vector<ShakeBond>& shake_bonds) {
    n_molecules_ = static_cast<int>(molecule_constraint_counts.size());
    n_constraints_ = static_cast<int>(shake_bonds.size());

    std::vector<int> offsets(n_molecules_ + 1, 0);
    for (int mol = 0; mol < n_molecules_; mol++) {
        offsets[mol + 1] = offsets[mol] + molecule_constraint_counts[mol];
    }
    if (offsets.back() != n_constraints_) {
        throw std::runtime_error("CUDA serial SHAKE molecule counts do not match the constraint list.");
    }

    molecule_constraint_offsets_ = HostDeviceBuffer<int>::from_vector(offsets, true);
    shake_bonds_ = HostDeviceBuffer<ShakeBond>::from_vector(shake_bonds, true);
    ready_ = std::make_unique<HostDeviceBuffer<std::uint8_t>>(n_constraints_, false, true);
    failed_ = std::make_unique<HostDeviceBuffer<int>>(1, true, true);
}

void CudaSerialConstraintSolver::apply(
    coord_t* d_coords,
    const coord_t* d_xcoords,
    const double* d_winv) {
    if (!enabled()) return;

    ready_->zero();
    failed_->zero();

    const int blocks = (n_molecules_ + kSerialShakeThreads - 1) / kSerialShakeThreads;
    serial_shake_by_molecule_kernel<<<blocks, kSerialShakeThreads>>>(
        n_molecules_,
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
    if (failed_->cpu_data_p[0]) std::exit(EXIT_FAILURE);
}

void CudaSerialShake::init_backend(Context& ctx) {
    const int* counts = data_.mol_n_shakes->cpu_data_p;
    std::vector<int> molecule_constraint_counts(counts, counts + ctx.n_molecules());
    const ShakeBond* bonds = data_.shake_bonds->cpu_data_p;
    std::vector<ShakeBond> shake_bonds(bonds, bonds + data_.n_constraints);
    solver_.init(molecule_constraint_counts, shake_bonds);
}

void CudaSerialShake::apply(Context& ctx, HostDeviceBuffer<coord_t>& xcoords) {
    apply_to(ctx, ctx.coords->gpu_data_p, xcoords.gpu_data_p);
}

void CudaSerialShake::apply_to(Context& ctx, coord_t* d_coords, coord_t* d_xcoords) {
    solver_.apply(d_coords, d_xcoords, ctx.winv->gpu_data_p);
}

void CudaSerialShake::initial_shake(Context& ctx) {
    HostDeviceBuffer<coord_t> xcoords(ctx.n_atoms, true, true);
    coord_t* d_coords = ctx.coords->gpu_data_p;
    coord_t* d_xcoords = xcoords.gpu_data_p;

    check_cuda(cudaMemcpy(
        d_xcoords,
        d_coords,
        ctx.n_atoms * sizeof(coord_t),
        cudaMemcpyDeviceToDevice));
    apply_to(ctx, d_coords, d_xcoords);

    ctx.coords->download();
    coord_t* coords = ctx.coords->cpu_data_p;
    vel_t* velocities = ctx.velocities->cpu_data_p;
    coord_t* host_xcoords = xcoords.cpu_data_p;
    for (int i = 0; i < ctx.n_atoms; i++) {
        host_xcoords[i].x = coords[i].x - ctx.dt * velocities[i].x;
        host_xcoords[i].y = coords[i].y - ctx.dt * velocities[i].y;
        host_xcoords[i].z = coords[i].z - ctx.dt * velocities[i].z;
    }

    xcoords.upload();
    apply_to(ctx, d_xcoords, d_coords);
    xcoords.download();

    for (int i = 0; i < ctx.n_atoms; i++) {
        velocities[i].x = (coords[i].x - host_xcoords[i].x) / ctx.dt;
        velocities[i].y = (coords[i].y - host_xcoords[i].y) / ctx.dt;
        velocities[i].z = (coords[i].z - host_xcoords[i].z) / ctx.dt;
    }
    ctx.velocities->upload();
}
