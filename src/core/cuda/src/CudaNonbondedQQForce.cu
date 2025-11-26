#include "cuda/include/CudaNonbondedQQForce.cuh"
#include "utils.h"
namespace CudaNonbondedQQForce {
bool is_initialized = false;
q_atom_t* d_q_atoms = nullptr;
q_charge_t* d_q_charges = nullptr;
int* d_LJ_matrix = nullptr;
bool* d_excluded = nullptr;
q_elscale_t* d_q_elscales = nullptr;
q_catype_t* d_q_catypes = nullptr;
q_atype_t* d_q_atypes = nullptr;
coord_t* d_coords = nullptr;
E_nonbonded_t* d_EQ_nonbond_qq = nullptr;
dvel_t* d_dvelocities = nullptr;
double* d_lambdas = nullptr;
}  // namespace CudaNonbondedQQForce

__global__ void calc_nonbonded_qq_forces_kernel(
    q_atom_t* q_atoms,
    q_charge_t* q_charges,
    int* LJ_matrix,
    bool* excluded,
    q_elscale_t* q_elscales,
    q_catype_t* q_catypes,
    q_atype_t* q_atypes,
    coord_t* coords,
    E_nonbonded_t* EQ_nonbond_qq,
    dvel_t* dvelocities,
    double* lambdas,
    int n_qatoms,
    int n_lambdas,
    int n_atoms_solute,
    topo_t topo,
    int n_qelscales

) {
    int idx = blockIdx.x * blockDim.x + threadIdx.x;
    if (idx >= n_qatoms) return;
    int qi = idx;

    int ai, aj;
    double crg_i, crg_j;
    double elscale, scaling;
    q_catype_t qi_type, qj_type;
    bool bond23, bond14;
    coord_t da;
    double r2a, ra, r6a;
    double Vela, V_a, V_b;
    double dva;
    double ai_aii, aj_aii, ai_bii, aj_bii;

    for (int state = 0; state < n_lambdas; state++) {
        for (int qj = qi + 1; qj < n_qatoms; qj++) {
            ai = q_atoms[qi].a - 1;
            aj = q_atoms[qj].a - 1;

            crg_i = q_charges[qi + n_qatoms * state].q;
            crg_j = q_charges[qj + n_qatoms * state].q;

            bond23 = LJ_matrix[ai * n_atoms_solute + aj] == 3;
            bond14 = LJ_matrix[ai * n_atoms_solute + aj] == 1;

            if (bond23) continue;
            if (excluded[ai] || excluded[aj]) continue;

            scaling = bond14 ? topo.el14_scale : 1;

            elscale = 1;
            for (int k = 0; k < n_qelscales; k++) {
                if (q_elscales[k + n_qelscales * state].qi == qi + 1 && q_elscales[k + n_qelscales * state].qj == qj + 1) {
                    elscale = q_elscales[k + n_qelscales * state].mu;
                }
            }

            qi_type = q_catypes[q_atypes[qi + n_qatoms * state].code - 1];
            qj_type = q_catypes[q_atypes[qj + n_qatoms * state].code - 1];

            da.x = coords[aj].x - coords[ai].x;
            da.y = coords[aj].y - coords[ai].y;
            da.z = coords[aj].z - coords[ai].z;
            r2a = 1 / (pow(da.x, 2) + pow(da.y, 2) + pow(da.z, 2));
            ra = sqrt(r2a);
            r6a = r2a * r2a * r2a;

            Vela = scaling * topo.coulomb_constant * crg_i * crg_j * ra * elscale;

            ai_aii = bond14 ? qi_type.Ai_14 : qi_type.Ai;
            aj_aii = bond14 ? qj_type.Ai_14 : qj_type.Ai;
            ai_bii = bond14 ? qi_type.Bi_14 : qi_type.Bi;
            aj_bii = bond14 ? qj_type.Bi_14 : qj_type.Bi;

            V_a = r6a * r6a * ai_aii * aj_aii;
            V_b = r6a * ai_bii * aj_bii;
            dva = r2a * (-Vela - 12 * V_a + 6 * V_b) * lambdas[state];

            atomicAdd(&dvelocities[ai].x, -dva * da.x);
            atomicAdd(&dvelocities[ai].y, -dva * da.y);
            atomicAdd(&dvelocities[ai].z, -dva * da.z);
            atomicAdd(&dvelocities[aj].x, dva * da.x);
            atomicAdd(&dvelocities[aj].y, dva * da.y);
            atomicAdd(&dvelocities[aj].z, dva * da.z);

            atomicAdd(&EQ_nonbond_qq[state].Ucoul, Vela);
            atomicAdd(&EQ_nonbond_qq[state].Uvdw, (V_a - V_b));
        }
    }
}

void calc_nonbonded_qq_forces_host() {
    using namespace CudaNonbondedQQForce;
    if (!is_initialized) {
        check_cudaMalloc((void**)&d_q_atoms, sizeof(q_atom_t) * n_qatoms);
        check_cudaMalloc((void**)&d_q_charges, sizeof(q_charge_t) * n_qatoms * n_lambdas);
        check_cudaMalloc((void**)&d_LJ_matrix, sizeof(int) * n_atoms_solute * n_atoms_solute);
        check_cudaMalloc((void**)&d_excluded, sizeof(bool) * n_atoms);
        check_cudaMalloc((void**)&d_q_elscales, sizeof(q_elscale_t) * n_qelscales);
        check_cudaMalloc((void**)&d_q_catypes, sizeof(q_catype_t) * n_qcatypes);
        check_cudaMalloc((void**)&d_q_atypes, sizeof(q_atype_t) * n_qatoms * n_lambdas);
        check_cudaMalloc((void**)&d_coords, sizeof(coord_t) * n_atoms);
        check_cudaMalloc((void**)&d_EQ_nonbond_qq, sizeof(E_nonbonded_t) * n_lambdas);
        check_cudaMalloc((void**)&d_dvelocities, sizeof(dvel_t) * n_atoms);
        check_cudaMalloc((void**)&d_lambdas, sizeof(double) * n_lambdas);

        cudaMemcpy(d_q_atoms, q_atoms, sizeof(q_atom_t) * n_qatoms, cudaMemcpyHostToDevice);
        cudaMemcpy(d_q_charges, q_charges, sizeof(q_charge_t) * n_qatoms * n_lambdas, cudaMemcpyHostToDevice);
        cudaMemcpy(d_LJ_matrix, LJ_matrix, sizeof(int) * n_atoms_solute * n_atoms_solute, cudaMemcpyHostToDevice);
        cudaMemcpy(d_excluded, excluded, sizeof(bool) * n_atoms, cudaMemcpyHostToDevice);
        cudaMemcpy(d_q_elscales, q_elscales, sizeof(q_elscale_t) * n_qelscales, cudaMemcpyHostToDevice);
        cudaMemcpy(d_q_catypes, q_catypes, sizeof(q_catype_t) * n_qcatypes, cudaMemcpyHostToDevice);
        cudaMemcpy(d_q_atypes, q_atypes, sizeof(q_atype_t) * n_qatoms * n_lambdas, cudaMemcpyHostToDevice);
        cudaMemcpy(d_lambdas, lambdas, sizeof(double) * n_lambdas, cudaMemcpyHostToDevice);

        is_initialized = true;
    }
    cudaMemcpy(d_coords, coords, sizeof(coord_t) * n_atoms, cudaMemcpyHostToDevice);
    cudaMemcpy(d_dvelocities, dvelocities, sizeof(dvel_t) * n_atoms, cudaMemcpyHostToDevice);
    cudaMemcpy(d_EQ_nonbond_qq, EQ_nonbond_qq, sizeof(E_nonbonded_t) * n_lambdas, cudaMemcpyHostToDevice);
    int blockSize = 256;
    int numBlocks = (n_qatoms + blockSize - 1) / blockSize;
    calc_nonbonded_qq_forces_kernel<<<numBlocks, blockSize>>>(
        d_q_atoms,
        d_q_charges,
        d_LJ_matrix,
        d_excluded,
        d_q_elscales,
        d_q_catypes,
        d_q_atypes,
        d_coords,
        d_EQ_nonbond_qq,
        d_dvelocities,
        d_lambdas,
        n_qatoms,
        n_lambdas,
        n_atoms_solute,
        topo, n_qelscales);
    cudaDeviceSynchronize();
    cudaMemcpy(EQ_nonbond_qq, d_EQ_nonbond_qq, sizeof(E_nonbonded_t) * n_lambdas, cudaMemcpyDeviceToHost);
    cudaMemcpy(dvelocities, d_dvelocities, sizeof(dvel_t) * n_atoms, cudaMemcpyDeviceToHost);
}
void cleanup_nonbonded_qq_force() {
    using namespace CudaNonbondedQQForce;
    if (is_initialized) {
        cudaFree(d_q_atoms);
        cudaFree(d_q_charges);
        cudaFree(d_LJ_matrix);
        cudaFree(d_excluded);
        cudaFree(d_q_elscales);
        cudaFree(d_q_catypes);
        cudaFree(d_q_atypes);
        cudaFree(d_coords);
        cudaFree(d_EQ_nonbond_qq);
        cudaFree(d_dvelocities);
        cudaFree(d_lambdas);
        is_initialized = false;
    }
}
