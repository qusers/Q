#include "constants.h"
#include "cuda_force_accumulation.cuh"
#include "cuda_water_boundary_force.cuh"
#include "geometry.h"

namespace {
const int kBlockSize = 256;

/*
 * Each thread processes one water molecule.
 *
 * Only the oxygen atom receives the radial force.
 */
__global__ void calc_radix_kernel(
    int n_waters,
    int n_atoms_solute,
    int n_atoms, int energy_stride,
    const coord_t* coords,
    const bool* excluded,
    const double* temperature_results,
    topo_t topo,
    md_t md,
    double Dwmz,
    double awmz,
    dvel_t* dvelocities,
    energy_accum_t* energy) {
    const int water_idx = blockIdx.x * blockDim.x + threadIdx.x;

    if (water_idx >= n_waters) return;

    const int oxygen_idx = n_atoms_solute + 3 * water_idx;

    if (excluded[oxygen_idx]) return;

    const int replica = blockIdx.y;
    const int atom_offset = replica * n_atoms;
    coords += atom_offset;
    temperature_results += replica * N_TRESULT;
    dvelocities += atom_offset;
    energy += replica * energy_stride;

    const double shift = md.radial_force != 0.0 ? sqrt(Boltz * temperature_results[R_TFREE] / md.radial_force) : 0.0;

    const double dx = coords[oxygen_idx].x - topo.solvent_center.x;

    const double dy = coords[oxygen_idx].y - topo.solvent_center.y;

    const double dz = coords[oxygen_idx].z - topo.solvent_center.z;

    const double b = sqrt(dx * dx + dy * dy + dz * dz);

    const double db = b - (topo.solvent_radius - shift);

    double ener;
    double dv;

    if (db > 0.0) {
        ener = 0.5 * md.radial_force * db * db - Dwmz;
        dv = md.radial_force * db / b;
    } else if (b > 0.0) {
        const double fexp = exp(awmz * db);
        ener = Dwmz * (fexp * fexp - 2.0 * fexp);
        dv = -2.0 * Dwmz * awmz * (fexp - fexp * fexp) / b;
    } else {
        ener = 0.0;
        dv = 0.0;
    }

    atomic_add_energy(&energy[E_RESTR_RADX], ener);

    atomic_add_force(&dvelocities[oxygen_idx].x, dv * dx);

    atomic_add_force(&dvelocities[oxygen_idx].y, dv * dy);

    atomic_add_force(&dvelocities[oxygen_idx].z, dv * dz);
}

__global__ void prepare_wshells_kernel(int n_shells, bool update_theta_corr, shell_t* wshells) {
    const int shell_idx = blockIdx.x * blockDim.x + threadIdx.x;
    if (shell_idx >= n_shells) return;

    const int replica = blockIdx.y;
    wshells += replica * n_shells;

    shell_t& shell = wshells[shell_idx];

    if (update_theta_corr) {
        shell.avtheta /= static_cast<double>(itdis_update);
        shell.avn_inshell /= static_cast<double>(itdis_update);
        shell.theta_corr += shell.avtheta - acos(shell.cstb);

        shell.avtheta = 0;
        shell.avn_inshell = 0;
    }

    shell.n_inshell = 0;
}

__global__ void calc_theta_and_shell_kernel(int n_waters,
                                            int n_atoms_solute,
                                            int n_shells,
                                            int n_max_inshell,
                                            const coord_t* coords,
                                            const bool* excluded,
                                            topo_t topo,
                                            shell_t* wshells,
                                            double* theta,
                                            int* list_sh) {
    const int water_idx = blockIdx.x * blockDim.x + threadIdx.x;

    if (water_idx >= n_waters) return;

    const int wi = n_atoms_solute + 3 * water_idx;

    if (excluded[wi]) return;

    const int replica = blockIdx.y;
    const int atom_offset = replica * (n_atoms_solute + 3 * n_waters);  // = replica * n_atoms
    coords += atom_offset;
    wshells += replica * n_shells;
    theta += replica * n_waters;
    list_sh += replica * n_shells * n_max_inshell;

    coord_t rmu = coords[wi + 1] + coords[wi + 2] - coords[wi] * 2.0;

    const double rm = norm(rmu);
    rmu = rmu / rm;

    coord_t rcu = coords[wi] - topo.solvent_center;

    const double rc = norm(rcu);
    rcu = rcu / rc;

    double cos_th = dot(rmu, rcu);

    if (cos_th > 1.0) {
        cos_th = 1.0;
    } else if (cos_th < -1.0) {
        cos_th = -1.0;
    }

    theta[water_idx] = acos(cos_th);

    const shell_t& inner_shell = wshells[n_shells - 1];

    const double inner_boundary = inner_shell.router - inner_shell.dr;

    if (rc <= inner_boundary) {
        return;
    }

    int shell_idx = n_shells - 1;

    while (shell_idx > 0 && rc > wshells[shell_idx].router) {
        --shell_idx;
    }

    const int position = atomicAdd(&wshells[shell_idx].n_inshell, 1);
    list_sh[shell_idx * n_max_inshell + position] = water_idx;
}

__global__ void calc_polx_force_with_rank_kernel(int n_atoms_solute,
                                                 int n_shells,
                                                 int n_max_inshell,
                                                 int n_atoms,
                                                 int n_waters,
                                                 int energy_stride,
                                                 const coord_t* coords,
                                                 const int* list_sh,
                                                 const double* theta,
                                                 shell_t* wshells,
                                                 topo_t topo,
                                                 md_t md,
                                                 dvel_t* dvelocities,
                                                 energy_accum_t* energy) {
    const int replica = blockIdx.y;
    const int atom_offset = replica * n_atoms;
    coords += atom_offset;
    list_sh += replica * n_shells * n_max_inshell;
    theta += replica * n_waters;
    wshells += replica * n_shells;
    dvelocities += atom_offset;
    energy += replica * energy_stride;

    const int flat_idx = blockIdx.x * blockDim.x + threadIdx.x;
    const int total_slots = n_shells * n_max_inshell;
    if (flat_idx >= total_slots) return;
    const int shell_idx = flat_idx / n_max_inshell;

    const int position = flat_idx % n_max_inshell;

    const int shell_size = wshells[shell_idx].n_inshell;

    if (position >= shell_size) return;

    const int shell_offset = shell_idx * n_max_inshell;

    const int water_idx = list_sh[shell_offset + position];

    if (water_idx < 0) return;

    const double current_theta = theta[water_idx];

    int rank = 0;
    for (int other_position = 0; other_position < shell_size; ++other_position) {
        const int other_water_idx = list_sh[shell_offset + other_position];

        const double other_theta = theta[other_water_idx];

        if (other_theta < current_theta || (other_theta == current_theta && other_water_idx < water_idx)) {
            ++rank;
        }
    }

    double cos_arg = 1.0 + (1.0 - 2.0 * (rank + 1)) / static_cast<double>(shell_size);

    if (cos_arg > 1.0) {
        cos_arg = 1.0;
    } else if (cos_arg < -1.0) {
        cos_arg = -1.0;
    }

    double target_theta = acos(cos_arg);

    target_theta -= 1.5 * sin(target_theta) * wshells[shell_idx].cstb;

    if (target_theta < 0.0) {
        target_theta = 0.0;
    } else if (target_theta > M_PI) {
        target_theta = M_PI;
    }

    const double dtheta = current_theta - target_theta + wshells[shell_idx].theta_corr;

    const double ener = 0.5 * md.polarisation_force * dtheta * dtheta;

    atomic_add_energy(&energy[E_RESTR_POLX], ener);

    const double dv = md.polarisation_force * dtheta;

    const int wi = n_atoms_solute + 3 * water_idx;

    coord_t rmu = coords[wi + 1] + coords[wi + 2] - coords[wi] * 2.0;

    const double rm = norm(rmu);
    rmu = rmu / rm;

    coord_t rcu = coords[wi] - topo.solvent_center;

    const double rc = norm(rcu);
    rcu = rcu / rc;

    double cos_th = dot(rmu, rcu);

    if (cos_th > 1.0) {
        cos_th = 1.0;
    } else if (cos_th < -1.0) {
        cos_th = -1.0;
    }

    double sin_th = sin(acos(cos_th));

    if (fabs(sin_th) < k_singular_sin_epsilon) {
        sin_th = k_singular_sin_epsilon;
    }

    const double f0 = -dv / sin_th;

    const coord_t g_mu = (rcu - rmu * cos_th) / rm;

    const coord_t g_c = (rmu - rcu * cos_th) / rc;

    const coord_t fO = g_mu * (-2.0) + g_c;

    const coord_t fH = g_mu;

    atomic_add_force(&dvelocities[wi].x, f0 * fO.x);

    atomic_add_force(&dvelocities[wi].y, f0 * fO.y);

    atomic_add_force(&dvelocities[wi].z, f0 * fO.z);

    atomic_add_force(&dvelocities[wi + 1].x, f0 * fH.x);

    atomic_add_force(&dvelocities[wi + 1].y, f0 * fH.y);

    atomic_add_force(&dvelocities[wi + 1].z, f0 * fH.z);

    atomic_add_force(&dvelocities[wi + 2].x, f0 * fH.x);

    atomic_add_force(&dvelocities[wi + 2].y, f0 * fH.y);

    atomic_add_force(&dvelocities[wi + 2].z, f0 * fH.z);

    /*
     * Equivalent to:
     *
     * avtheta += sum(theta) / shell_size
     */
    atomicAdd(&wshells[shell_idx].avtheta, current_theta / static_cast<double>(shell_size));

    /*
     * One thread per non-empty shell updates avn_inshell.
     */
    if (position == 0) {
        atomicAdd(&wshells[shell_idx].avn_inshell, static_cast<double>(shell_size));
    }
}

}  // namespace

void CudaWaterBoundaryForce::calc(Context& ctx, int iteration) {
    if (ctx.n_waters() <= 0) return;
    calc_radix(ctx);

    if (ctx.md.polarisation) {
        calc_polx(ctx, iteration);
    }
}

void CudaWaterBoundaryForce::calc_radix(
    Context& ctx) {
    const int num_blocks =
        (ctx.n_waters() + kBlockSize - 1) /
        kBlockSize;

    calc_radix_kernel<<<dim3(num_blocks, ctx.n_replicates()), kBlockSize>>>(
        ctx.n_waters(),
        ctx.n_atoms_solute, ctx.n_atoms, ctx.energy.replica_stride(),
        ctx.coords->gpu_data_p,
        ctx.excluded->gpu_data_p,
        temperature_->data().results->gpu_data_p,
        ctx.topo,
        ctx.md,
        data_.Dwmz,
        data_.awmz,
        ctx.dvelocities->gpu_data_p,
        ctx.energy.device());
}

void CudaWaterBoundaryForce::calc_polx(Context& ctx, int iteration) {
    const int n_shells = data_.n_shells;
    const int n_max_inshell = data_.n_max_inshell;

    if (n_shells <= 0 || n_max_inshell <= 0 || ctx.n_waters() <= 0) return;

    shell_t* device_wshells = data_.wshells->gpu_data_p;
    double* device_theta = data_.theta->gpu_data_p;
    int* device_list_sh = data_.list_sh->gpu_data_p;

    const bool update_theta_corr = iteration != 0 && iteration % itdis_update == 0;
    const int shell_blocks = (n_shells + kBlockSize - 1) / kBlockSize;

    prepare_wshells_kernel<<<dim3(shell_blocks, ctx.n_replicates()), kBlockSize>>>(n_shells, update_theta_corr, device_wshells);

    const int water_blocks = (ctx.n_waters() + kBlockSize - 1) / kBlockSize;
    calc_theta_and_shell_kernel<<<dim3(water_blocks, ctx.n_replicates()), kBlockSize>>>(ctx.n_waters(),
                                                                                        ctx.n_atoms_solute,
                                                                                        n_shells,
                                                                                        n_max_inshell,
                                                                                        ctx.coords->gpu_data_p,
                                                                                        ctx.excluded->gpu_data_p,
                                                                                        ctx.topo,
                                                                                        device_wshells,
                                                                                        device_theta,
                                                                                        device_list_sh);
    check_cuda(cudaGetLastError());

    const int total_slots = n_shells * n_max_inshell;
    const int force_blocks = (total_slots + kBlockSize - 1) / kBlockSize;
    calc_polx_force_with_rank_kernel<<<dim3(force_blocks, ctx.n_replicates()), kBlockSize>>>(
        ctx.n_atoms_solute,
        n_shells,
        n_max_inshell,
        ctx.n_atoms,
        ctx.n_waters(),
        ctx.energy.replica_stride(),
        ctx.coords->gpu_data_p,
        device_list_sh,
        device_theta,
        device_wshells,
        ctx.topo,
        ctx.md,
        ctx.dvelocities->gpu_data_p,
        ctx.energy.device());
}

void CudaWaterBoundaryForce::sync_for_output(
    Context& ctx) {
    if (data_.wshells &&
        ctx.md.polarisation) {
        data_.wshells->download();
    }

    WaterBoundaryForce::sync_for_output(ctx);
}
