#include "constants.h"
#include "common/include/context.h"
#include "common/include/coulomb.h"
#include "common/include/softcore.h"
#include "common/include/vdw_rules.h"
#include "cuda/include/cuda_force_accum.cuh"
#include "cuda/include/cuda_nonbonded_14_force.cuh"
#include "cuda_utility.cuh"

namespace CudaNonbonded14Force {
bool is_initialized = false;
constexpr int kNonbonded14ModeCount = 3;

int* d_atom_to_qi = nullptr;
real_t* d_evdw_totals = nullptr;
real_t* d_ecoul_totals = nullptr;

__device__ __forceinline__ nonbond_work_t nonbond14_rsqrt(nonbond_work_t value) {
#ifdef QDYN_SPFP
    return rsqrtf(value);
#else
    return rsqrt(value);
#endif
}

__device__ __forceinline__ int unified_parameter_index(
    int atom_idx,
    int state,
    int n_atoms,
    int n_qatoms,
    const int* atom_to_qi) {
    const int qi = atom_to_qi[atom_idx];
    if (qi == -1 || state == 0) {
        return atom_idx;
    }
    return n_atoms + (state - 1) * n_qatoms + qi;
}
// todo: lambda only has 2, so we can simplify it and only call once.

__device__ void calculate_nonbonded_14_pair(
    const coord_t& x,
    const coord_t& y,
    real_t x_charge,
    real_t y_charge,
    real_t x_aii_normal,
    real_t y_aii_normal,
    real_t x_bii_normal,
    real_t y_bii_normal,
    real_t x_aii,
    real_t y_aii,
    real_t x_bii,
    real_t y_bii,
    nonbond_work_t coulomb_constant,
    nonbond_work_t scaling,
    int vdw_rule,
    nonbond_work_t lambda,
    real_t softcore_alpha,
    bool softcore_use_max_potential,
    nonbond_work_t& evdw,
    nonbond_work_t& ecoul,
    nonbond_work_t& dv) {
    const nonbond_work_t dx = static_cast<nonbond_work_t>(x.x - y.x);
    const nonbond_work_t dy = static_cast<nonbond_work_t>(x.y - y.y);
    const nonbond_work_t dz = static_cast<nonbond_work_t>(x.z - y.z);
    const nonbond_work_t r = nonbond14_rsqrt(dx * dx + dy * dy + dz * dz);
    const nonbond_work_t r2 = r * r;
    const nonbond_work_t r6 = r2 * r2 * r2;

    ecoul = fortran_coulomb_energy(
        static_cast<nonbond_work_t>(x_charge),
        static_cast<nonbond_work_t>(y_charge),
        r,
        coulomb_constant,
        scaling) * lambda;

    nonbond_work_t pair_a = 0.0;
    nonbond_work_t pair_b = 0.0;
    if (vdw_rule == VDW_GEOMETRIC) {
        calc_vdw_geometric(
            static_cast<nonbond_work_t>(x_aii),
            static_cast<nonbond_work_t>(y_aii),
            static_cast<nonbond_work_t>(x_bii),
            static_cast<nonbond_work_t>(y_bii),
            static_cast<nonbond_work_t>(1.0),
            &pair_a,
            &pair_b);
    } else {
        calc_vdw_arithmetic(
            static_cast<nonbond_work_t>(x_aii),
            static_cast<nonbond_work_t>(y_aii),
            static_cast<nonbond_work_t>(x_bii),
            static_cast<nonbond_work_t>(y_bii),
            static_cast<nonbond_work_t>(1.0),
            &pair_a,
            &pair_b);
    }
    nonbond_work_t lookup_pair_a = pair_a;
    nonbond_work_t lookup_pair_b = pair_b;
    if (softcore_use_max_potential) {
        if (vdw_rule == VDW_GEOMETRIC) {
            calc_vdw_geometric(
                static_cast<nonbond_work_t>(x_aii_normal),
                static_cast<nonbond_work_t>(y_aii_normal),
                static_cast<nonbond_work_t>(x_bii_normal),
                static_cast<nonbond_work_t>(y_bii_normal),
                static_cast<nonbond_work_t>(1.0),
                &lookup_pair_a,
                &lookup_pair_b);
        } else {
            calc_vdw_arithmetic(
                static_cast<nonbond_work_t>(x_aii_normal),
                static_cast<nonbond_work_t>(y_aii_normal),
                static_cast<nonbond_work_t>(x_bii_normal),
                static_cast<nonbond_work_t>(y_bii_normal),
                static_cast<nonbond_work_t>(1.0),
                &lookup_pair_a,
                &lookup_pair_b);
        }
    }
    const nonbond_work_t softcore_lookup = q_softcore_lookup_value(
        static_cast<nonbond_work_t>(softcore_alpha),
        lookup_pair_a,
        lookup_pair_b,
        softcore_use_max_potential);
    nonbond_work_t v_a = 0.0;
    nonbond_work_t v_b = 0.0;
    nonbond_work_t force_lj_term = 0.0;
    q_softcore_lj(pair_a, pair_b, r6, softcore_lookup, &v_a, &v_b, &force_lj_term);
    evdw = (v_a - v_b) * lambda;
    dv = r2 * (-ecoul - force_lj_term * lambda);
}

__global__ void calc_nonbonded_14_force_kernel(
    int n_pairs,
    const int3* ngbrs_14,
    topo_t d_topo,
    const bool* d_excluded,
    const int* atom_to_qi,
    const ccharge_t* unified_ccharges,
    const catype_t* unified_catypes,
    const coord_t* d_coords,
    cuda_dvel_t* d_dvelocities,
    real_t* evdw_totals,
    real_t* ecoul_totals,
    bool include_pp,
    bool include_qq,
    int state,
    int n_atoms,
    int n_qatoms,
    real_t lambda,
    bool softcore_enabled,
    const real_t* q_softcore_values,
    bool softcore_use_max_potential) {
    const int idx = blockIdx.x * blockDim.x + threadIdx.x;
    if (idx >= n_pairs) return;

    const int3 pair = ngbrs_14[idx];
    const int ai = pair.x;
    const int aj = pair.y;
    const int mode = pair.z;

    if (ai < 0 || aj < 0 || ai >= n_atoms || aj >= n_atoms || ai == aj) return;
    if (d_excluded[ai] || d_excluded[aj]) return;
    if (mode < NONBONDED_14_PP || mode > NONBONDED_14_QQ) return;
    if (include_pp) {
        if (mode != NONBONDED_14_PP) return;
    } else if (mode == NONBONDED_14_PP) {
        return;
    }
    if (!include_qq && mode == NONBONDED_14_QQ) return;

    const int ai_param_idx = unified_parameter_index(ai, state, n_atoms, n_qatoms, atom_to_qi);
    const int aj_param_idx = unified_parameter_index(aj, state, n_atoms, n_qatoms, atom_to_qi);

    const ccharge_t ai_charge = unified_ccharges[ai_param_idx];
    const ccharge_t aj_charge = unified_ccharges[aj_param_idx];
    const catype_t ai_type = unified_catypes[ai_param_idx];
    const catype_t aj_type = unified_catypes[aj_param_idx];
    const coord_t ri = d_coords[ai];
    const coord_t rj = d_coords[aj];

    nonbond_work_t evdw = 0.0;
    nonbond_work_t ecoul = 0.0;
    nonbond_work_t dv = 0.0;
    const nonbond_work_t pair_lambda = static_cast<nonbond_work_t>((mode == NONBONDED_14_PP) ? 1.0 : lambda);
    real_t softcore_alpha = static_cast<real_t>(0.0);
    if (softcore_enabled && mode != NONBONDED_14_PP) {
        const int ai_qi = atom_to_qi[ai];
        const int aj_qi = atom_to_qi[aj];
        const real_t ai_alpha = (ai_qi >= 0) ? q_softcore_values[ai_qi] : static_cast<real_t>(0.0);
        const real_t aj_alpha = (aj_qi >= 0) ? q_softcore_values[aj_qi] : static_cast<real_t>(0.0);
        if (ai_qi >= 0 && aj_qi >= 0) {
            softcore_alpha = q_softcore_pair_value(ai_alpha, aj_alpha, softcore_use_max_potential);
        } else if (ai_qi >= 0) {
            softcore_alpha = ai_alpha;
        } else if (aj_qi >= 0) {
            softcore_alpha = aj_alpha;
        }
    }

    calculate_nonbonded_14_pair(
        ri,
        rj,
        ai_charge.charge,
        aj_charge.charge,
        ai_type.aii_normal,
        aj_type.aii_normal,
        ai_type.bii_normal,
        aj_type.bii_normal,
        ai_type.aii_1_4,
        aj_type.aii_1_4,
        ai_type.bii_1_4,
        aj_type.bii_1_4,
        static_cast<nonbond_work_t>(d_topo.coulomb_constant),
        static_cast<nonbond_work_t>(d_topo.el14_scale),
        d_topo.vdw_rule,
        pair_lambda,
        softcore_alpha,
        softcore_use_max_potential,
        evdw,
        ecoul,
        dv);

    const nonbond_work_t dx = static_cast<nonbond_work_t>(rj.x - ri.x);
    const nonbond_work_t dy = static_cast<nonbond_work_t>(rj.y - ri.y);
    const nonbond_work_t dz = static_cast<nonbond_work_t>(rj.z - ri.z);
    atomic_add_force_component(&d_dvelocities[ai].x, static_cast<real_t>(-dv * dx));
    atomic_add_force_component(&d_dvelocities[ai].y, static_cast<real_t>(-dv * dy));
    atomic_add_force_component(&d_dvelocities[ai].z, static_cast<real_t>(-dv * dz));
    atomic_add_force_component(&d_dvelocities[aj].x, static_cast<real_t>(dv * dx));
    atomic_add_force_component(&d_dvelocities[aj].y, static_cast<real_t>(dv * dy));
    atomic_add_force_component(&d_dvelocities[aj].z, static_cast<real_t>(dv * dz));

    atomicAdd(&evdw_totals[mode], evdw);
    atomicAdd(&ecoul_totals[mode], ecoul);
}

}  // namespace CudaNonbonded14Force

namespace {
struct Nonbonded14EnergyBuckets {
    real_t evdw[CudaNonbonded14Force::kNonbonded14ModeCount] = {};
    real_t ecoul[CudaNonbonded14Force::kNonbonded14ModeCount] = {};
};
}

static Nonbonded14EnergyBuckets calc_nonbonded_14_force_state_host(
    int state,
    real_t lambda,
    bool include_pp,
    bool include_qq) {
    using namespace CudaNonbonded14Force;

    auto& host = Context::instance();
    const int n_ngbrs_14 = host.n_ngbrs14;
    Nonbonded14EnergyBuckets energies = {};
    if (n_ngbrs_14 == 0) return energies;

    cudaMemset(d_ecoul_totals, 0, sizeof(real_t) * kNonbonded14ModeCount);
    cudaMemset(d_evdw_totals, 0, sizeof(real_t) * kNonbonded14ModeCount);

    const int block_size = 256;
    const int num_blocks = (n_ngbrs_14 + block_size - 1) / block_size;
    const bool softcore_enabled = host.n_qsoftcores > 0 && !include_pp;
    const real_t* q_softcore_values = softcore_enabled
        ? host.q_softcore_values->gpu_data_p + state * host.n_qatoms
        : nullptr;

    calc_nonbonded_14_force_kernel<<<num_blocks, block_size>>>(
        n_ngbrs_14,
        host.ngbrs_14->gpu_data_p,
        host.topo,
        host.excluded->gpu_data_p,
        d_atom_to_qi,
        host.unified_ccharges->gpu_data_p,
        host.unified_catypes->gpu_data_p,
        host.coords->gpu_data_p,
        cuda_force_accum_buffer(host),
        d_evdw_totals,
        d_ecoul_totals,
        include_pp,
        include_qq,
        state,
        host.n_atoms,
        host.n_qatoms,
        lambda,
        softcore_enabled,
        q_softcore_values,
        host.q_softcore_use_max_potential);

    cudaDeviceSynchronize();

    cudaMemcpy(energies.evdw, d_evdw_totals, sizeof(real_t) * kNonbonded14ModeCount, cudaMemcpyDeviceToHost);
    cudaMemcpy(energies.ecoul, d_ecoul_totals, sizeof(real_t) * kNonbonded14ModeCount, cudaMemcpyDeviceToHost);

    return energies;
}

void calc_nonbonded_14_forces_host() {
    auto& host = Context::instance();

    if (host.n_ngbrs14 == 0) return;

    Nonbonded14EnergyBuckets pp_energies = calc_nonbonded_14_force_state_host(0, 1.0, true, false);
    host.E_nonbond_pp.Uvdw += pp_energies.evdw[NONBONDED_14_PP];
    host.E_nonbond_pp.Ucoul += pp_energies.ecoul[NONBONDED_14_PP];

    if (host.n_lambdas == 0) return;

    auto *lambdas = host.lambdas->cpu_data_p;
    for (int state = 0; state < host.n_lambdas; state++) {
        const real_t lambda = lambdas[state];
        // CUDA Q-Q nonbonds, including Q-Q 1-4 pairs, are handled by the CPU
        // fallback in CudaHandler to keep this correctness path aligned with Fortran.
        Nonbonded14EnergyBuckets energies = calc_nonbonded_14_force_state_host(state, lambda, false, false);

        if (lambda != 0.0) {
            host.EQ_nonbond_qp[state].Uvdw += energies.evdw[NONBONDED_14_QP] / lambda;
            host.EQ_nonbond_qp[state].Ucoul += energies.ecoul[NONBONDED_14_QP] / lambda;
            host.EQ_nonbond_qq[state].Uvdw += energies.evdw[NONBONDED_14_QQ] / lambda;
            host.EQ_nonbond_qq[state].Ucoul += energies.ecoul[NONBONDED_14_QQ] / lambda;
        }
    }
}

void init_nonbonded_14_force_kernel_data() {
    using namespace CudaNonbonded14Force;
    if (is_initialized) return;

    const auto& host = Context::instance();

    check_cudaMalloc((void**)&d_atom_to_qi, sizeof(int) * host.atom_to_qi.size());
    check_cuda(cudaMemcpy(d_atom_to_qi, host.atom_to_qi.data(), sizeof(int) * host.atom_to_qi.size(), cudaMemcpyHostToDevice));

    check_cudaMalloc((void**)&d_evdw_totals, sizeof(real_t) * kNonbonded14ModeCount);
    check_cudaMalloc((void**)&d_ecoul_totals, sizeof(real_t) * kNonbonded14ModeCount);

    is_initialized = true;
}

void cleanup_nonbonded_14_force() {
    using namespace CudaNonbonded14Force;
    if (!is_initialized) return;

    cudaFree(d_atom_to_qi);
    cudaFree(d_evdw_totals);
    cudaFree(d_ecoul_totals);

    d_atom_to_qi = nullptr;
    d_evdw_totals = nullptr;
    d_ecoul_totals = nullptr;
    is_initialized = false;
}
