#include "constants.h"
#include "common/include/context.h"
#include "common/include/vdw_rules.h"
#include "cuda/include/cuda_context.cuh"
#include "cuda/include/cuda_nonbonded_14_force.cuh"
#include "cuda_utility.cuh"

namespace CudaNonbonded14Force {
bool is_initialized = false;
constexpr int kNonbonded14ModeCount = 3;

int* d_atom_to_qi = nullptr;
double* d_evdw_totals = nullptr;
double* d_ecoul_totals = nullptr;

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
    double x_charge,
    double y_charge,
    double x_aii,
    double y_aii,
    double x_bii,
    double y_bii,
    double coulomb_constant,
    double scaling,
    int vdw_rule,
    double lambda,
    double& evdw,
    double& ecoul,
    double& dv) {
    const double3 d = {x.x - y.x, x.y - y.y, x.z - y.z};
    const double r = rsqrt(d.x * d.x + d.y * d.y + d.z * d.z);
    const double r2 = r * r;
    const double r6 = r2 * r2 * r2;

    ecoul = scaling * coulomb_constant * x_charge * y_charge * r * lambda;

    double v_a = 0.0;
    double v_b = 0.0;
    if (vdw_rule == VDW_GEOMETRIC) {
        calc_vdw_geometric(x_aii, y_aii, x_bii, y_bii, r6, &v_a, &v_b);
    } else {
        calc_vdw_arithmetic(x_aii, y_aii, x_bii, y_bii, r6, &v_a, &v_b);
    }
    v_a *= lambda;
    v_b *= lambda;
    evdw = v_a - v_b;
    dv = r2 * (-ecoul - 12.0 * v_a + 6.0 * v_b);
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
    dvel_t* d_dvelocities,
    double* evdw_totals,
    double* ecoul_totals,
    bool include_pp,
    int state,
    int n_atoms,
    int n_qatoms,
    double lambda) {
    const int idx = blockIdx.x * blockDim.x + threadIdx.x;
    if (idx >= n_pairs) return;

    const int3 pair = ngbrs_14[idx];
    const int ai = pair.x;
    const int aj = pair.y;
    const int mode = pair.z;

    if (ai < 0 || aj < 0 || ai >= n_atoms || aj >= n_atoms || ai == aj) return;
    if (d_excluded[ai] || d_excluded[aj]) return;
    if (mode < NONBONDED_14_PP || mode > NONBONDED_14_QQ) return;
    if (!include_pp && mode == NONBONDED_14_PP) return;

    const int ai_param_idx = unified_parameter_index(ai, state, n_atoms, n_qatoms, atom_to_qi);
    const int aj_param_idx = unified_parameter_index(aj, state, n_atoms, n_qatoms, atom_to_qi);

    const ccharge_t ai_charge = unified_ccharges[ai_param_idx];
    const ccharge_t aj_charge = unified_ccharges[aj_param_idx];
    const catype_t ai_type = unified_catypes[ai_param_idx];
    const catype_t aj_type = unified_catypes[aj_param_idx];
    const coord_t ri = d_coords[ai];
    const coord_t rj = d_coords[aj];

    double evdw = 0.0;
    double ecoul = 0.0;
    double dv = 0.0;
    const double pair_lambda = (mode == NONBONDED_14_PP) ? 1.0 : lambda;

    calculate_nonbonded_14_pair(
        ri,
        rj,
        ai_charge.charge,
        aj_charge.charge,
        ai_type.aii_1_4,
        aj_type.aii_1_4,
        ai_type.bii_1_4,
        aj_type.bii_1_4,
        d_topo.coulomb_constant,
        d_topo.el14_scale,
        d_topo.vdw_rule,
        pair_lambda,
        evdw,
        ecoul,
        dv);

    const double3 d = {rj.x - ri.x, rj.y - ri.y, rj.z - ri.z};
    atomicAdd(&d_dvelocities[ai].x, -dv * d.x);
    atomicAdd(&d_dvelocities[ai].y, -dv * d.y);
    atomicAdd(&d_dvelocities[ai].z, -dv * d.z);
    atomicAdd(&d_dvelocities[aj].x, dv * d.x);
    atomicAdd(&d_dvelocities[aj].y, dv * d.y);
    atomicAdd(&d_dvelocities[aj].z, dv * d.z);

    atomicAdd(&evdw_totals[mode], evdw);
    atomicAdd(&ecoul_totals[mode], ecoul);
}

}  // namespace CudaNonbonded14Force

namespace {
struct Nonbonded14EnergyBuckets {
    double evdw[CudaNonbonded14Force::kNonbonded14ModeCount] = {};
    double ecoul[CudaNonbonded14Force::kNonbonded14ModeCount] = {};
};
}

static Nonbonded14EnergyBuckets calc_nonbonded_14_force_state_host(
    int state,
    double lambda,
    bool include_pp) {
    using namespace CudaNonbonded14Force;

    auto& host = Context::instance();
    auto& context = CudaContext::instance();
    const int n_ngbrs_14 = host.n_ngbrs14;
    Nonbonded14EnergyBuckets energies = {};
    if (n_ngbrs_14 == 0) return energies;

    cudaMemset(d_ecoul_totals, 0, sizeof(double) * kNonbonded14ModeCount);
    cudaMemset(d_evdw_totals, 0, sizeof(double) * kNonbonded14ModeCount);

    const int block_size = 256;
    const int num_blocks = (n_ngbrs_14 + block_size - 1) / block_size;

    calc_nonbonded_14_force_kernel<<<num_blocks, block_size>>>(
        n_ngbrs_14,
        context.d_ngbrs_14,
        host.topo,
        context.d_excluded,
        d_atom_to_qi,
        context.d_unified_ccharges,
        context.d_unified_catypes,
        context.d_coords,
        context.d_dvelocities,
        d_evdw_totals,
        d_ecoul_totals,
        include_pp,
        state,
        host.n_atoms,
        host.n_qatoms,
        lambda);

    cudaDeviceSynchronize();

    cudaMemcpy(energies.evdw, d_evdw_totals, sizeof(double) * kNonbonded14ModeCount, cudaMemcpyDeviceToHost);
    cudaMemcpy(energies.ecoul, d_ecoul_totals, sizeof(double) * kNonbonded14ModeCount, cudaMemcpyDeviceToHost);

    return energies;
}

void calc_nonbonded_14_forces_host() {
    auto& host = Context::instance();

    if (host.n_ngbrs14 == 0) return;

    for (int state = 0; state < host.n_lambdas; state++) {
        const double lambda = host.lambdas[state];
        const bool include_pp = (state == 0);
        Nonbonded14EnergyBuckets energies = calc_nonbonded_14_force_state_host(state, lambda, include_pp);

        if (include_pp) {
            host.E_nonbond_pp.Uvdw += energies.evdw[NONBONDED_14_PP];
            host.E_nonbond_pp.Ucoul += energies.ecoul[NONBONDED_14_PP];
        }

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

    check_cudaMalloc((void**)&d_evdw_totals, sizeof(double) * kNonbonded14ModeCount);
    check_cudaMalloc((void**)&d_ecoul_totals, sizeof(double) * kNonbonded14ModeCount);

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
