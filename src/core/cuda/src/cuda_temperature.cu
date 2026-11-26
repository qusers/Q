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

__global__ void calc_temperature_kernel(int n_atoms, int n_atoms_solute, const double* masses, vel_t* velocities, bool* excluded, double boltz, double ekinmax,
                                        energy_accum_t* Temp_solute, energy_accum_t* Tfree_solute, energy_accum_t* Texcl_solute,
                                        energy_accum_t* Temp_solvent, energy_accum_t* Tfree_solvent, energy_accum_t* Texcl_solvent) {
    const int idx = blockIdx.x * blockDim.x + threadIdx.x;

    // Kinetic energy of this atom in fixed point, in one of four bins:
    // free solute, excluded solute, free solvent, excluded solvent. Threads past the last atom add 0.
    energy_accum_t bin[4] = {0, 0, 0, 0};
    if (idx < n_atoms) {
        double mass_i = masses[idx];
        const double vx = velocities[idx].x;
        const double vy = velocities[idx].y;
        const double vz = velocities[idx].z;
        double ener = .5 * mass_i * (vx * vx + vy * vy + vz * vz);

        bin[(idx < n_atoms_solute ? 0 : 2) + (excluded[idx] ? 1 : 0)] = energy_to_accum(ener);
        if (ener > ekinmax) {
            printf(">>> WARNING: hot atom %d: %f\n", idx, ener / boltz / 3);
        }
    }

    // Sum over the warp, then one atomic add per slot and warp instead of two per atom. The sums are
    // integers, so the totals do not depend on the order.
    const unsigned mask = 0xffffffffu;
    for (int offset = 16; offset > 0; offset >>= 1) {
        for (int b = 0; b < 4; ++b) {
            bin[b] += __shfl_down_sync(mask, bin[b], offset);
        }
    }
    if ((threadIdx.x & 31) != 0) return;

    auto add = [](energy_accum_t* target, energy_accum_t value) {
        if (value != 0) atomicAdd(reinterpret_cast<unsigned long long int*>(target), static_cast<unsigned long long int>(value));
    };
    add(Temp_solute, bin[0] + bin[1]);
    add(Tfree_solute, bin[0]);
    add(Texcl_solute, bin[1]);
    add(Temp_solvent, bin[2] + bin[3]);
    add(Tfree_solvent, bin[2]);
    add(Texcl_solvent, bin[3]);
}

__global__ void finalize_temperature_kernel(
    const energy_accum_t* acc,
    int Ndegf, int Ndegfree, int Ndegfree_solute, int Ndegfree_solvent,
    bool separate_scaling, double tau_T, double dt, double target_T,
    double* results) {
    if (threadIdx.x || blockIdx.x) return;

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
    double Uexcl = energy_from_accum(acc[TEXCL_SOL]) + energy_from_accum(acc[TEXCL_SLV]);
    double ts_excl = 1.0;
    if (Ndegf > Ndegfree) {
        double Texcl = 2.0 * Uexcl / Boltz / (Ndegf - Ndegfree);
        if (Texcl > 0) ts_excl = sqrt(1 + (dt / tau_T) * (target_T / Texcl - 1.0));
    }
    results[R_TSCALE_EXCL] = ts_excl;

    results[R_TEMP] = Temp;
    results[R_TFREE] = Tfree;
    results[R_TSCALE_SOL] = ts_sol;
    results[R_TSCALE_SLV] = ts_slv;
    results[R_UKIN] = Ukin;
}

}  // namespace

void CudaTemperature::init_backend(Context& ctx) {
    accum_ = std::make_unique<HostDeviceBuffer<energy_accum_t>>(N_TACCUM);
}

void CudaTemperature::calc(Context& ctx) {
    accum_->zero();  // resets host + device slots
    energy_accum_t* d = accum_->gpu_data_p;

    double Ekinmax = 1000.0 * data_.Ndegf * Boltz * ctx.md.temperature / 2.0 / ctx.n_atoms;
    int blockSize = 256;
    int numBlocks = (ctx.n_atoms + blockSize - 1) / blockSize;

    calc_temperature_kernel<<<numBlocks, blockSize>>>(
        ctx.n_atoms, ctx.n_atoms_solute,
        ctx.masses->gpu_data_p,
        ctx.velocities->gpu_data_p, ctx.excluded->gpu_data_p, Boltz, Ekinmax,
        d + TEMP_SOL, d + TFREE_SOL, d + TEXCL_SOL,
        d + TEMP_SLV, d + TFREE_SLV, d + TEXCL_SLV);

    finalize_temperature_kernel<<<1, 1>>>(
        d, data_.Ndegf, data_.Ndegfree, data_.Ndegfree_solute, data_.Ndegfree_solvent,
        data_.separate_scaling, data_.tau_T, ctx.dt, ctx.md.temperature,
        data_.results->gpu_data_p);
}

void CudaTemperature::sync_for_output(Context& ctx) {
    data_.results->download();  // one small memcpy on log steps
    Temperature::sync_for_output(ctx);
}