#include "constants.h"
#include "cuda_force_accumulation.cuh"
#include "cuda_temperature.cuh"
#include "md_types.h"

namespace {
enum TAccum { TEMP_SOL,
              TFREE_SOL,
              TEXCL_SOL,
              TEMP_SLV,
              TFREE_SLV,
              TEXCL_SLV,
              N_TACCUM };

__global__ void calc_temperature_kernel(int n_atoms, int n_atoms_solute, atype_t* atypes, catype_t* catypes, vel_t* velocities, bool* excluded, double boltz, double ekinmax,
                                        energy_accum_t* Temp_solute, energy_accum_t* Tfree_solute, energy_accum_t* Texcl_solute,
                                        energy_accum_t* Temp_solvent, energy_accum_t* Tfree_solvent, energy_accum_t* Texcl_solvent) {
    int idx = blockIdx.x * blockDim.x + threadIdx.x;
    const int replica = blockIdx.y;
    if (idx >= n_atoms) return;

    velocities += replica * n_atoms;

    const int accum_offset = replica * N_TACCUM;
    Temp_solute += accum_offset;
    Tfree_solute += accum_offset;
    Texcl_solute += accum_offset;
    Temp_solvent += accum_offset;
    Tfree_solvent += accum_offset;
    Texcl_solvent += accum_offset;

    double mass_i = catypes[atypes[idx].code - 1].m;
    const double vx = velocities[idx].x;
    const double vy = velocities[idx].y;
    const double vz = velocities[idx].z;
    double ener = .5 * mass_i * (vx * vx + vy * vy + vz * vz);

    bool is_solute = (idx < n_atoms_solute);
    bool is_excluded = excluded[idx];
    if (is_solute) {
        atomic_add_energy(Temp_solute, ener);
        if (!is_excluded)
            atomic_add_energy(Tfree_solute, ener);
        else
            atomic_add_energy(Texcl_solute, ener);
    } else {
        atomic_add_energy(Temp_solvent, ener);
        if (!is_excluded)
            atomic_add_energy(Tfree_solvent, ener);
        else
            atomic_add_energy(Texcl_solvent, ener);
    }
    if (ener > ekinmax) {
        printf(">>> WARNING: hot atom %d: %f\n", idx, ener / boltz / 3);
    }
}

__global__ void finalize_temperature_kernel(
    const energy_accum_t* acc,
    int Ndegf, int Ndegfree, int Ndegfree_solute, int Ndegfree_solvent,
    bool separate_scaling, double tau_T, double dt, double target_T,
    double* results) {
    if (threadIdx.x != 0) return;
    const int replica = blockIdx.x;
    results += replica * N_TRESULT;
    acc += replica * N_TACCUM;

    double Tf_sol = energy_from_accum(acc[TFREE_SOL]);
    double Tf_slv = energy_from_accum(acc[TFREE_SLV]);
    double Ukin = energy_from_accum(acc[TEMP_SOL]) + energy_from_accum(acc[TEMP_SLV]);

    double Temp = 2.0 * Ukin / Boltz / Ndegf;
    double Tfree = 2.0 * (Tf_sol + Tf_slv) / Boltz / Ndegfree;

    double ts_sol, ts_slv;
    if (separate_scaling) {
        double t_slv = 2.0 * Tf_slv / Boltz / Ndegfree_solvent;
        double t_sol = 2.0 * Tf_sol / Boltz / Ndegfree_solute;
        ts_slv = (t_slv != 0) ? sqrt(1 + (dt / tau_T) * (target_T / t_slv - 1.0)) : 1.0;
        ts_sol = (t_sol != 0) ? sqrt(1 + (dt / tau_T) * (target_T / t_sol - 1.0)) : 1.0;
    } else {
        ts_slv = (Tfree != 0) ? sqrt(1 + (dt / tau_T) * (target_T / Tfree - 1.0)) : 1.0;
        ts_sol = ts_slv;
    }

    results[R_TEMP] = Temp;
    results[R_TFREE] = Tfree;
    results[R_TSCALE_SOL] = ts_sol;
    results[R_TSCALE_SLV] = ts_slv;
    results[R_UKIN] = Ukin;
}

}  // namespace

void CudaTemperature::init_backend(Context& ctx) {
    accum_ = std::make_unique<HostDeviceBuffer<energy_accum_t>>(N_TACCUM * ctx.n_replicates());
}

void CudaTemperature::calc(Context& ctx) {
    accum_->zero();  // resets host + device slots
    energy_accum_t* d = accum_->gpu_data_p;

    double Ekinmax = 1000.0 * data_.Ndegf * Boltz * ctx.md.temperature / 2.0 / ctx.n_atoms;
    int blockSize = 256;
    int numBlocks = (ctx.n_atoms + blockSize - 1) / blockSize;

    calc_temperature_kernel<<<dim3(numBlocks, ctx.n_replicates()), blockSize>>>(
        ctx.n_atoms, ctx.n_atoms_solute,
        ctx.atypes->gpu_data_p, ctx.catypes->gpu_data_p,
        ctx.velocities->gpu_data_p, ctx.excluded->gpu_data_p, Boltz, Ekinmax,
        d + TEMP_SOL, d + TFREE_SOL, d + TEXCL_SOL,
        d + TEMP_SLV, d + TFREE_SLV, d + TEXCL_SLV);

    finalize_temperature_kernel<<<ctx.n_replicates(), 1>>>(
        d, data_.Ndegf, data_.Ndegfree, data_.Ndegfree_solute, data_.Ndegfree_solvent,
        data_.separate_scaling, data_.tau_T, ctx.dt, ctx.md.temperature,
        data_.results->gpu_data_p);
}

void CudaTemperature::sync_for_output(Context& ctx, int replica) {
    if (replica == 0) {
        data_.results->download();  // one small memcpy on log steps
    }
    Temperature::sync_for_output(ctx, replica);
}