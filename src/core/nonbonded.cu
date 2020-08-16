#include "system.h"
#include "nonbonded.h"
#include "utils.h"
#include <stdio.h>

/* =============================================
 * == NON-BONDED INTERACTIONS
 * =============================================
 */

coord_t *X;
dvel_t *DV_X;
dvel_t *PP_MAT, *h_PP_MAT;
bool pp_gpu_set = false;

// Constants pointers
ccharge_t *D_ccharges;
charge_t *D_charges;
catype_t *D_catypes;
atype_t * D_atypes;
p_atom_t *D_patoms;
int *D_LJ_matrix;
bool *D_excluded;

void calc_nonbonded_pp_forces() {
    bool bond14, bond23;
    double scaling;
    coord_t da;
    double r2a, ra, r6a;
    double Vela, V_a, V_b;
    double dva;
    double crg_i, crg_j;
    double ai_aii, aj_aii, ai_bii, aj_bii;
    int i, j;
    catype_t ai_type, aj_type;

    for (int pi = 0; pi < n_patoms; pi++) {
        for (int pj = pi+1; pj < n_patoms; pj++) {
            i = p_atoms[pi].a - 1;
            j = p_atoms[pj].a - 1;
            bond23 = LJ_matrix[i * n_atoms_solute + j] == 3;
            bond14 = LJ_matrix[i * n_atoms_solute + j] == 1;

            if (bond23) continue;
            if (excluded[i] || excluded[j]) continue;

            scaling = bond14 ? .5 : 1;

            crg_i = ccharges[charges[i].code - 1].charge;
            crg_j = ccharges[charges[j].code - 1].charge;

            ai_type = catypes[atypes[i].code - 1];
            aj_type = catypes[atypes[j].code - 1];

            da.x = coords[j].x - coords[i].x;
            da.y = coords[j].y - coords[i].y;
            da.z = coords[j].z - coords[i].z;
            r2a = 1 / (pow(da.x, 2) + pow(da.y, 2) + pow(da.z, 2));
            ra = sqrt(r2a);
            r6a = r2a * r2a * r2a;

            Vela = scaling * Coul * crg_i * crg_j * ra;

            ai_aii = bond14 ? ai_type.aii_1_4 : ai_type.aii_normal;
            aj_aii = bond14 ? aj_type.aii_1_4 : aj_type.aii_normal;
            ai_bii = bond14 ? ai_type.bii_1_4 : ai_type.bii_normal;
            aj_bii = bond14 ? aj_type.bii_1_4 : aj_type.bii_normal;

            V_a = r6a * r6a * ai_aii * aj_aii;
            V_b = r6a * ai_bii * aj_bii;
            dva = r2a * ( -Vela -12 * V_a + 6 * V_b);

            dvelocities[i].x -= dva * da.x;
            dvelocities[i].y -= dva * da.y;
            dvelocities[i].z -= dva * da.z;

            dvelocities[j].x += dva * da.x;
            dvelocities[j].y += dva * da.y;
            dvelocities[j].z += dva * da.z;

            energies.Ucoul += Vela;
            energies.Uvdw += (V_a - V_b);
        }
    }
}

void calc_nonbonded_pp_forces_host() {
    int mem_size_X = n_atoms_solute * sizeof(coord_t);
    int mem_size_DV_X = n_atoms_solute * sizeof(dvel_t);

    int mem_size_ccharges = n_ccharges * sizeof(ccharge_t);
    int mem_size_charges = n_atoms * sizeof(charge_t);
    int mem_size_catypes = n_catypes * sizeof(catype_t);
    int mem_size_atypes = n_atoms * sizeof(atype_t);
    int mem_size_patoms = n_patoms * sizeof(p_atom_t);
    int mem_size_LJ_matrix = n_atoms_solute * n_atoms_solute * sizeof(int);
    int mem_size_excluded = n_atoms * sizeof(bool);

    int mem_size_PP_MAT = n_patoms * n_patoms * sizeof(dvel_t);

    if (!pp_gpu_set) {
        pp_gpu_set = true;
        #ifdef DEBUG
        printf("Allocating P\n");
        #endif
        check_cudaMalloc((void**) &X, mem_size_X);
        #ifdef DEBUG
        printf("Allocating DV_P\n");
        #endif
        check_cudaMalloc((void**) &DV_X, mem_size_DV_X);

        #ifdef DEBUG
        printf("Allocating D_ccharges\n");
        #endif
        check_cudaMalloc((void**) &D_ccharges, mem_size_ccharges);
        #ifdef DEBUG
        printf("Allocating D_charges\n");
        #endif
        check_cudaMalloc((void**) &D_charges, mem_size_charges);
        #ifdef DEBUG
        printf("Allocating D_catypes\n");
        #endif
        check_cudaMalloc((void**) &D_catypes, mem_size_catypes);
        #ifdef DEBUG
        printf("Allocating D_atypes\n");
        #endif
        check_cudaMalloc((void**) &D_atypes, mem_size_atypes);
        #ifdef DEBUG
        printf("Allocating D_patoms\n");
        #endif
        check_cudaMalloc((void**) &D_patoms, mem_size_patoms);
        #ifdef DEBUG
        printf("Allocating D_LJ_matrix\n");
        #endif
        check_cudaMalloc((void**) &D_LJ_matrix, mem_size_LJ_matrix);
        #ifdef DEBUG
        printf("Allocating D_excluded\n");
        #endif
        check_cudaMalloc((void**) &D_excluded, mem_size_excluded);

        #ifdef DEBUG
        printf("Allocating PP_MAT\n");
        #endif
        check_cudaMalloc((void**) &PP_MAT, mem_size_PP_MAT);

        cudaMemcpy(D_ccharges, ccharges, mem_size_ccharges, cudaMemcpyHostToDevice);
        cudaMemcpy(D_charges, charges, mem_size_charges, cudaMemcpyHostToDevice);
        cudaMemcpy(D_catypes, catypes, mem_size_catypes, cudaMemcpyHostToDevice);
        cudaMemcpy(D_atypes, atypes, mem_size_atypes, cudaMemcpyHostToDevice);
        cudaMemcpy(D_patoms, p_atoms, mem_size_patoms, cudaMemcpyHostToDevice);
        cudaMemcpy(D_LJ_matrix, LJ_matrix, mem_size_LJ_matrix, cudaMemcpyHostToDevice);
        cudaMemcpy(D_excluded, excluded, mem_size_excluded, cudaMemcpyHostToDevice);

        // for (int i = 0; i < 100; i++) {
        //     for (int j = 0; j < 100; j++) {
        //         printf("LJ_matrix[%d][%d] = %d\n", i, j, LJ_matrix[i * n_patoms + j]);
        //     }
        // }

        h_PP_MAT = (dvel_t*) malloc(mem_size_PP_MAT);
    }

    cudaMemcpy(X, coords, mem_size_X, cudaMemcpyHostToDevice);
    cudaMemcpy(DV_X, dvelocities, mem_size_DV_X, cudaMemcpyHostToDevice);

    dim3 threads,grid;

    threads = dim3(BLOCK_SIZE, BLOCK_SIZE);
    grid = dim3((n_patoms + BLOCK_SIZE - 1) / threads.x, (n_patoms + BLOCK_SIZE - 1) / threads.y);
    double evdw, ecoul;

    calc_pp_dvel_matrix<<<grid, threads>>>(n_patoms, n_atoms_solute, X, &evdw, &ecoul, PP_MAT, D_ccharges, D_charges, D_catypes, D_atypes, D_patoms, D_LJ_matrix, D_excluded);
    calc_pp_dvel_vector<<<((n_patoms+BLOCK_SIZE - 1) / BLOCK_SIZE), BLOCK_SIZE>>>(n_patoms, DV_X, PP_MAT, D_patoms);

    cudaMemcpy(dvelocities, DV_X, mem_size_DV_X, cudaMemcpyDeviceToHost);

    #ifdef DEBUG
    cudaMemcpy(h_PP_MAT, PP_MAT, mem_size_PP_MAT, cudaMemcpyDeviceToHost);
    #endif

    // for (int i = 0; i < n_patoms; i++) {
    //     printf("excluded[%d] = %s\n", i, excluded[i] ? "True" : "False");
    // }

    #ifdef DEBUG
    for (int i = 0; i < n_patoms; i++) {
        for (int j = 0; j < n_patoms; j++) {
            if (h_PP_MAT[i * n_patoms + j].x > 100)
            printf("PP_MAT[%d][%d] = %f %f %f\n", i, j, h_PP_MAT[i * n_patoms + j].x, h_PP_MAT[i * n_patoms + j].y, h_PP_MAT[i * n_patoms + j].z);
        }
    }
    #endif
}

/* =============================================
 * == DEVICE
 * =============================================
 */

__device__ void calc_pp_dvel_matrix_incr(int row, int pi, int column, int pj,
    coord_t *Xs, coord_t *Ys, int *LJs, bool *excluded_s, double *Evdw, double *Ecoul, dvel_t *patom_a, dvel_t *patom_b,
    ccharge_t *D_ccharges, charge_t *D_charges, catype_t *D_catypes, atype_t *D_atypes, p_atom_t *D_patoms) {
    
    bool bond14, bond23;
    double scaling;
    coord_t da;
    double r2a, ra, r6a;
    double Vela, V_a, V_b;
    double dva;
    double crg_i, crg_j;
    double ai_aii, aj_aii, ai_bii, aj_bii;
    catype_t ai_type, aj_type;

    bond23 = LJs[row * BLOCK_SIZE + column] == 3;
    bond14 = LJs[row * BLOCK_SIZE + column] == 1;

    if (bond23) return;
    if (excluded_s[row] || excluded_s[BLOCK_SIZE + column]) return;

    scaling = bond14 ? .5 : 1;

    crg_i = D_ccharges[D_charges[pi].code - 1].charge;
    crg_j = D_ccharges[D_charges[pj].code - 1].charge;

    ai_type = D_catypes[D_atypes[pi].code - 1];
    aj_type = D_catypes[D_atypes[pj].code - 1];

    da.x = Ys[column].x - Xs[row].x;
    da.y = Ys[column].y - Xs[row].y;
    da.z = Ys[column].z - Xs[row].z;
    r2a = 1 / (pow(da.x, 2) + pow(da.y, 2) + pow(da.z, 2));
    ra = sqrt(r2a);
    r6a = r2a * r2a * r2a;

    Vela = scaling * Coul * crg_i * crg_j * ra;

    ai_aii = bond14 ? ai_type.aii_1_4 : ai_type.aii_normal;
    aj_aii = bond14 ? aj_type.aii_1_4 : aj_type.aii_normal;
    ai_bii = bond14 ? ai_type.bii_1_4 : ai_type.bii_normal;
    aj_bii = bond14 ? aj_type.bii_1_4 : aj_type.bii_normal;

    V_a = r6a * r6a * ai_aii * aj_aii;
    V_b = r6a * ai_bii * aj_bii;
    dva = r2a * ( -Vela -12 * V_a + 6 * V_b);

    patom_a->x = -dva * da.x;
    patom_a->y = -dva * da.y;
    patom_a->z = -dva * da.z;

    patom_b->x = dva * da.x;
    patom_b->y = dva * da.y;
    patom_b->z = dva * da.z;

    *Evdw += Vela;
    *Ecoul += (V_a - V_b);
}

__global__ void calc_pp_dvel_matrix(int n_patoms, int n_atoms_solute,
    coord_t *X, double *Evdw, double *Ecoul, dvel_t *PP_MAT,
    ccharge_t *D_ccharges, charge_t *D_charges, catype_t *D_catypes, atype_t *D_atypes, p_atom_t *D_patoms, int *D_LJ_matrix, bool *D_excluded) {
    // Block index
    int bx = blockIdx.x;
    int by = blockIdx.y;

    // Thread index
    int tx = threadIdx.x;
    int ty = threadIdx.y;

    // if (bx == 9 && by == 0) printf("bx = %d by = %d tx = %d ty = %d\n", bx, by, tx, ty);

    int aStart = BLOCK_SIZE * by;
    int bStart = BLOCK_SIZE * bx;

    if (aStart + ty >= n_patoms) return;
    if (bStart + tx >= n_patoms) return;

    __shared__ coord_t Xs[BLOCK_SIZE];
    __shared__ coord_t Ys[BLOCK_SIZE];
    __shared__     int LJs[BLOCK_SIZE * BLOCK_SIZE];
    __shared__    bool excluded_s[2 * BLOCK_SIZE];
    
    // if (bx == 0 && by == 0) printf("bx = %d by = %d tx = %d ty = %d\n", bx, by, tx, ty);

    int pi = D_patoms[aStart + ty].a-1;
    int pj = D_patoms[bStart + tx].a-1;

    if (tx == 0) {
        Xs[ty] = X[pi];
        excluded_s[ty] = D_excluded[pi];
    }

    if (ty == 0) {
        Ys[tx] = X[pj];
        excluded_s[BLOCK_SIZE + tx] = D_excluded[pj];
    }
    LJs[ty * BLOCK_SIZE + tx] = D_LJ_matrix[pi * n_atoms_solute + pj];

    // if (bx == 0 && by == 1 && tx == 0) printf("excluded_s[%d] = %s\n", ty, excluded_s[ty] ? "True" : "False");
    // if (bx == 0 && by == 1 && ty == 0) printf("excluded_s[%d] = %s\n", BLOCK_SIZE + tx, excluded_s[BLOCK_SIZE + tx] ? "True" : "False");

    __syncthreads();

    if (bx < by || (bx == by && tx < ty)) return;

    dvel_t patom_a, patom_b;
    memset(&patom_a, 0, sizeof(dvel_t));
    memset(&patom_b, 0, sizeof(dvel_t));

    int row = by * BLOCK_SIZE + ty;
    int column = bx * BLOCK_SIZE + tx;

    if (bx != by || tx != ty) {
        double evdw, ecoul;
        // __device__ void calc_pp_dvel_matrix_incr(int row, int pi, int column, int pj,
        //     coord_t *Xs, coord_t *Ys, bool *LJs, double *Evdw, double *Ecoul, dvel_t *patom_a, dvel_t *patom_b,
        //     ccharge_t *D_ccharges, charge_t *D_charges, catype_t *D_catypes, atype_t *D_atypes, p_atom_t *D_patoms);

        calc_pp_dvel_matrix_incr(ty, pi, tx, pj, Xs, Ys, LJs, excluded_s, &evdw, &ecoul, &patom_a,
             &patom_b, D_ccharges, D_charges, D_catypes, D_atypes, D_patoms);
    }

    // if (row == 0 && column == 1) {
    //     printf("water_a = %f %f %f water_b = %f %f %f\n", water_a[0].x, water_a[0].y, water_a[0].z, water_b[0].x, water_b[0].y, water_b[0].z);
    // }

    PP_MAT[row * n_patoms + column] = patom_a;
    PP_MAT[column * n_patoms + row] = patom_b;

    __syncthreads();
}

__global__ void calc_pp_dvel_vector(int n_patoms, dvel_t *DV_X, dvel_t *PP_MAT, p_atom_t *D_patoms) {
    int row = blockIdx.x*blockDim.x + threadIdx.x;
    if (row >= n_patoms) return;

    dvel_t dP;
    dP.x = 0;
    dP.y = 0;
    dP.z = 0;

    for (int i = 0; i < n_patoms; i++) {
        if (i != row) {
            dP.x += PP_MAT[i + n_patoms * row].x;
            dP.y += PP_MAT[i + n_patoms * row].y;
            dP.z += PP_MAT[i + n_patoms * row].z;
        }
    }

    int p = D_patoms[row].a-1;

    DV_X[p].x += dP.x;
    DV_X[p].y += dP.y;
    DV_X[p].z += dP.z;

    __syncthreads();
}

void clean_d_patoms() {
    cudaFree(X);
    cudaFree(DV_X);
    cudaFree(PP_MAT);

    free(h_PP_MAT);

    cudaFree(D_ccharges);
    cudaFree(D_charges);
    cudaFree(D_catypes);
    cudaFree(D_atypes);
    cudaFree(D_patoms);
}