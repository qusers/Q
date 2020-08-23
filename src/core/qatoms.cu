// TODO: Add impropers, bond pairs
#include "system.h"
#include "qatoms.h"
#include "utils.h"
#include <math.h>
#include <stdio.h>

// Device pointers
calc_qw_t *QW_MAT, *h_QW_MAT;
calc_qp_t *QP_MAT, *h_QP_MAT;
dvel_t *QQ_MAT;

double *D_QP_Evdw, *D_QP_Ecoul, *h_QP_Evdw, *h_QP_Ecoul;
double *D_QP_evdw_TOT, *D_QP_ecoul_TOT, QP_evdw_TOT, QP_ecoul_TOT;

double *D_QW_Evdw, *D_QW_Ecoul, *h_QW_Evdw, *h_QW_Ecoul;
double *D_QW_evdw_TOT, *D_QW_ecoul_TOT, QW_evdw_TOT, QW_ecoul_TOT;

// Constants pointers
q_catype_t *D_qcatypes;
q_atype_t *D_qatypes;
q_charge_t *D_qcharges;
q_atom_t *D_qatoms;
double *D_lambdas;

bool qp_gpu_set = false;

void calc_nonbonded_qp_forces() {
    int i, j;
    coord_t da;
    double r2, r6, r;
    double ai_aii, aj_aii, ai_bii, aj_bii;
    q_catype_t qi_type;
    catype_t aj_type;
    bool bond23, bond14;
    double scaling, Vel, V_a, V_b, dv;

    for (int qi = 0; qi < n_qatoms; qi++) {
        for (int pj = 0; pj < n_patoms; pj++) {
            i = q_atoms[qi].a - 1;
            j = p_atoms[pj].a - 1;

            bond23 = LJ_matrix[i * n_atoms_solute + j] == 3;
            bond14 = LJ_matrix[i * n_atoms_solute + j] == 1;

            if (bond23) continue;
            if (excluded[i] || excluded[j]) continue;

            scaling = bond14 ? topo.el14_scale : 1;

            da.x = coords[j].x - coords[i].x;
            da.y = coords[j].y - coords[i].y;
            da.z = coords[j].z - coords[i].z;

            r2 = pow(da.x, 2) + pow(da.y, 2) + pow(da.z, 2);

            r6 = r2 * r2 * r2;
            r2 = 1 / r2;
            r = sqrt(r2);

            for (int state = 0; state < n_lambdas; state++) {
                qi_type = q_catypes[q_atypes[qi + n_qatoms * state].code - 1];
                aj_type = catypes[atypes[j].code - 1];

                ai_aii = bond14 ? qi_type.Ai_14 : qi_type.Ai;
                aj_aii = bond14 ? aj_type.aii_1_4 : aj_type.aii_normal;
                ai_bii = bond14 ? qi_type.Bi_14 : qi_type.Bi;
                aj_bii = bond14 ? aj_type.bii_1_4 : aj_type.bii_normal;

                Vel = topo.coulomb_constant * scaling * q_charges[qi + n_qatoms * state].q * ccharges[charges[j].code - 1].charge * r;
                V_a = ai_aii * aj_aii / (r6 * r6);
                V_b = ai_bii * aj_bii / r6;
                dv = r2 * (-Vel - (12 * V_a - 6 * V_b)) * lambdas[state];

                // Update forces
                dvelocities[i].x -= dv * da.x;
                dvelocities[i].y -= dv * da.y;
                dvelocities[i].z -= dv * da.z;
                dvelocities[j].x += dv * da.x;
                dvelocities[j].y += dv * da.y;
                dvelocities[j].z += dv * da.z;

                // Update Q totals
                EQ_nonbond_qp[state].Ucoul += Vel;
                EQ_nonbond_qp[state].Uvdw += (V_a - V_b);
            }
        }
    }
}

void calc_nonbonded_qp_forces_host() {
    int mem_size_X = n_atoms_solute * sizeof(coord_t);
    int mem_size_DV_X = n_atoms_solute * sizeof(dvel_t);
    int mem_size_QP_MAT = n_qatoms * n_patoms * sizeof(calc_qp_t);
    int mem_size_QQ_MAT = n_qatoms * n_qatoms * sizeof(dvel_t);

    int mem_size_qcatypes = n_qcatypes * sizeof(q_catype_t);
    int mem_size_qatypes = n_qatoms * n_lambdas * sizeof(q_atype_t);
    int mem_size_qcharges = n_qatoms * n_lambdas * sizeof(q_charge_t);
    int mem_size_qatoms = n_qatoms * sizeof(q_atom_t);
    int mem_size_lambdas = n_lambdas * sizeof(double);

    int mem_size_ccharges = n_ccharges * sizeof(ccharge_t);
    int mem_size_charges = n_atoms * sizeof(charge_t);
    int mem_size_catypes = n_catypes * sizeof(catype_t);
    int mem_size_atypes = n_atoms * sizeof(atype_t);
    int mem_size_patoms = n_patoms * sizeof(p_atom_t);
    int mem_size_LJ_matrix = n_atoms_solute * n_atoms_solute * sizeof(int);
    int mem_size_excluded = n_atoms * sizeof(bool);

    int n_blocks_q = (n_qatoms + BLOCK_SIZE - 1) / BLOCK_SIZE;
    int n_blocks_p = (n_patoms + BLOCK_SIZE - 1) / BLOCK_SIZE;
    //TODO make Evdw & Ecoul work for # of states > 2
    int mem_size_QP_Evdw = min(n_lambdas, 2) * n_blocks_q * n_blocks_p * sizeof(double);
    int mem_size_QP_Ecoul = min(n_lambdas, 2) * n_blocks_q * n_blocks_p * sizeof(double);

    if (!qp_gpu_set){
        qp_gpu_set = true;
        #ifdef DEBUG
        printf("Allocating X\n");
        #endif
        check_cudaMalloc((void**) &X, mem_size_X);
        #ifdef DEBUG
        printf("Allocating DV_X\n");
        #endif
        check_cudaMalloc((void**) &DV_X, mem_size_DV_X);

        #ifdef DEBUG
        printf("Allocating D_qcatypes\n");
        #endif
        check_cudaMalloc((void**) &D_qcatypes, mem_size_qcatypes);
        #ifdef DEBUG
        printf("Allocating D_qatypes\n");
        #endif
        check_cudaMalloc((void**) &D_qatypes, mem_size_qatypes);
        #ifdef DEBUG
        printf("Allocating D_qcharges\n");
        #endif
        check_cudaMalloc((void**) &D_qcharges, mem_size_qcharges);
        #ifdef DEBUG
        printf("Allocating D_qatoms\n");
        #endif
        check_cudaMalloc((void**) &D_qatoms, mem_size_qatoms);
        #ifdef DEBUG
        printf("Allocating D_lambdas\n");
        #endif
        check_cudaMalloc((void**) &D_lambdas, mem_size_lambdas);

        cudaMemcpy(D_qcatypes, q_catypes, mem_size_qcatypes, cudaMemcpyHostToDevice);
        cudaMemcpy(D_qatypes, q_atypes, mem_size_qatypes, cudaMemcpyHostToDevice);
        cudaMemcpy(D_qcharges, q_charges, mem_size_qcharges, cudaMemcpyHostToDevice);
        cudaMemcpy(D_qatoms, q_atoms, mem_size_qatoms, cudaMemcpyHostToDevice);
        cudaMemcpy(D_lambdas, lambdas, mem_size_lambdas, cudaMemcpyHostToDevice);

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

        cudaMemcpy(D_ccharges, ccharges, mem_size_ccharges, cudaMemcpyHostToDevice);
        cudaMemcpy(D_charges, charges, mem_size_charges, cudaMemcpyHostToDevice);
        cudaMemcpy(D_catypes, catypes, mem_size_catypes, cudaMemcpyHostToDevice);
        cudaMemcpy(D_atypes, atypes, mem_size_atypes, cudaMemcpyHostToDevice);
        cudaMemcpy(D_patoms, p_atoms, mem_size_patoms, cudaMemcpyHostToDevice);
        cudaMemcpy(D_LJ_matrix, LJ_matrix, mem_size_LJ_matrix, cudaMemcpyHostToDevice);
        cudaMemcpy(D_excluded, excluded, mem_size_excluded, cudaMemcpyHostToDevice);

        #ifdef DEBUG
        printf("Allocating QP_MAT\n");
        #endif
        check_cudaMalloc((void**) &QP_MAT, mem_size_QP_MAT);
        #ifdef DEBUG
        printf("Allocating QP_MAT\n");
        #endif
        check_cudaMalloc((void**) &QQ_MAT, mem_size_QQ_MAT);

        #ifdef DEBUG
        printf("Allocating D_QP_Evdw\n");
        #endif
        check_cudaMalloc((void**) &D_QP_Evdw, mem_size_QP_Evdw);    
        #ifdef DEBUG
        printf("Allocating D_QP_Ecoul\n");
        #endif
        check_cudaMalloc((void**) &D_QP_Ecoul, mem_size_QP_Ecoul);

        check_cudaMalloc((void**) &D_QP_evdw_TOT, sizeof(double));
        check_cudaMalloc((void**) &D_QP_ecoul_TOT, sizeof(double)); 

        h_QP_Evdw = (double*) malloc(mem_size_QP_Evdw);
        h_QP_Ecoul = (double*) malloc(mem_size_QP_Ecoul);

        h_QP_MAT = (calc_qp_t*) malloc(mem_size_QP_MAT);
    }

    cudaMemcpy(X, coords, mem_size_X, cudaMemcpyHostToDevice);
    cudaMemcpy(DV_X, dvelocities, mem_size_DV_X, cudaMemcpyHostToDevice);

    dim3 threads,grid;

    threads = dim3(BLOCK_SIZE, BLOCK_SIZE);
    grid = dim3((n_patoms + BLOCK_SIZE - 1) / threads.x, (n_qatoms + BLOCK_SIZE - 1) / threads.y);

    calc_qp_dvel_matrix<<<grid, threads>>>(n_qatoms, n_patoms, n_lambdas, n_atoms_solute, X, D_QP_Evdw, D_QP_Ecoul, QP_MAT,
        D_qcatypes, D_qatypes, D_qcharges, D_patoms, D_qatoms, D_lambdas, D_LJ_matrix, D_excluded,
        D_catypes, D_atypes, D_ccharges, D_charges, topo);
    calc_qp_dvel_vector_column<<<((n_patoms+BLOCK_SIZE - 1) / BLOCK_SIZE), BLOCK_SIZE>>>(n_qatoms, n_patoms, DV_X, QP_MAT, D_patoms);
    calc_qp_dvel_vector_row<<<((n_qatoms+BLOCK_SIZE - 1) / BLOCK_SIZE), BLOCK_SIZE>>>(n_qatoms, n_patoms, DV_X, QP_MAT, D_qatoms);

    #ifdef DEBUG
    cudaMemcpy(h_QP_MAT, QP_MAT, mem_size_QP_MAT, cudaMemcpyDeviceToHost);
    #endif

    #ifdef DEBUG
    for (int i = 0; i < n_qatoms; i++) {
        for (int j = 0; j < n_patoms; j++) {
            if (i == 0)
            // if (h_QP_MAT[i * n_patoms + j].Q.x > 100)
            printf("QP_MAT[%d][%d].Q = %f %f %f\n", i, j, h_QP_MAT[i * n_patoms + j].Q.x, h_QP_MAT[i * n_patoms + j].Q.y, h_QP_MAT[i * n_patoms + j].Q.z);
        }
    }
    #endif

    cudaMemcpy(dvelocities, DV_X, mem_size_DV_X, cudaMemcpyDeviceToHost);

    //TODO make Evdw & Ecoul work for # of states > 2
    for (int state = 0; state < min(2, n_lambdas); state++) {
        calc_energy_sum<<<1, threads>>>(n_blocks_q, n_blocks_p, D_QP_evdw_TOT, D_QP_ecoul_TOT, &D_QP_Evdw[state * n_blocks_p * n_blocks_q], &D_QP_Ecoul[state * n_blocks_p * n_blocks_q], false);

        cudaMemcpy(&QP_evdw_TOT, D_QP_evdw_TOT, sizeof(double), cudaMemcpyDeviceToHost);
        cudaMemcpy(&QP_ecoul_TOT, D_QP_ecoul_TOT, sizeof(double), cudaMemcpyDeviceToHost);
    
        EQ_nonbond_qp[state].Uvdw += QP_evdw_TOT;
        EQ_nonbond_qp[state].Ucoul += QP_ecoul_TOT;
    }
}

void calc_nonbonded_qw_forces() {
    int i;
    coord_t dO, dH1, dH2;
    double r2O, rH1, rH2, r6O, rO, r2H1, r2H2;
    double dvO, dvH1, dvH2;
    double V_a, V_b, VelO, VelH1, VelH2;
    q_catype_t qi_type;
    double ai_aii, ai_bii;

    if (A_O == 0) {
        catype_t catype_ow;    // Atom type of first O, H atom

        catype_ow = catypes[atypes[n_atoms_solute].code - 1];

        A_O = catype_ow.aii_normal;
        B_O = catype_ow.bii_normal;
    }

    // Loop over O-atoms, q-atoms
    for (int j = n_atoms_solute; j < n_atoms; j+= 3) {
        for (int qi = 0; qi < n_qatoms; qi++) {
            i = q_atoms[qi].a - 1;
            if (excluded[i] || excluded[j]) continue;
            dO.x = coords[j].x - coords[i].x;
            dO.y = coords[j].y - coords[i].y;
            dO.z = coords[j].z - coords[i].z;
            dH1.x = coords[j+1].x - coords[i].x;
            dH1.y = coords[j+1].y - coords[i].y;
            dH1.z = coords[j+1].z - coords[i].z;
            dH2.x = coords[j+2].x - coords[i].x;
            dH2.y = coords[j+2].y - coords[i].y;
            dH2.z = coords[j+2].z - coords[i].z;
            r2O = pow(dO.x, 2) + pow(dO.y, 2) + pow(dO.z, 2);
            rH1 = sqrt(1.0 / (pow(dH1.x, 2) + pow(dH1.y, 2) + pow(dH1.z, 2)));
            rH2 = sqrt(1.0 / (pow(dH2.x, 2) + pow(dH2.y, 2) + pow(dH2.z, 2)));
            r6O = r2O * r2O * r2O;
            r2O = 1.0 / r2O;
            rO = sqrt(r2O);
            r2H1 = rH1 * rH1;
            r2H2 = rH2 * rH2;

            // Reset potential
            dvO = 0;
            dvH1 = 0;
            dvH2 = 0;

            for (int state = 0; state < n_lambdas; state++) {
                qi_type = q_catypes[q_atypes[qi + n_qatoms * state].code - 1];

                ai_aii = qi_type.Ai;
                ai_bii = qi_type.Bi;

                V_a = ai_aii * A_O / (r6O * r6O);
                V_b = ai_bii * B_O / (r6O);

                VelO = topo.coulomb_constant * crg_ow * q_charges[qi + n_qatoms * state].q * rO;
                VelH1 = topo.coulomb_constant * crg_hw * q_charges[qi + n_qatoms * state].q * rH1;
                VelH2 = topo.coulomb_constant * crg_hw * q_charges[qi + n_qatoms * state].q * rH2;

                // if (state == 0 && qi == 1) printf("j = %d ai__aii = %f A_O = %f B_O = %f V_a = %f V_b = %f r6O = %f\n", j, ai_aii, A_O, B_O, V_a, V_b, r6O);

                dvO += r2O * (-VelO - (12 * V_a - 6 * V_b)) * lambdas[state];
                dvH1 -= r2H1 * VelH1 * lambdas[state];
                dvH2 -= r2H2 * VelH2 * lambdas[state];

                EQ_nonbond_qw[state].Ucoul += (VelO + VelH1 + VelH2);
                EQ_nonbond_qw[state].Uvdw += (V_a - V_b);
            }

            // Note r6O is not the usual 1/rO^6, but rather rO^6. be careful!!!

            // Update forces on Q-atom
            dvelocities[i].x -= (dvO * dO.x + dvH1 * dH1.x + dvH2 * dH2.x);
            dvelocities[i].y -= (dvO * dO.y + dvH1 * dH1.y + dvH2 * dH2.y);
            dvelocities[i].z -= (dvO * dO.z + dvH1 * dH1.z + dvH2 * dH2.z);

            // Update forces on water
            dvelocities[j].x += dvO * dO.x;
            dvelocities[j].y += dvO * dO.y;
            dvelocities[j].z += dvO * dO.z;
            dvelocities[j+1].x += dvH1 * dH1.x;
            dvelocities[j+1].y += dvH1 * dH1.y;
            dvelocities[j+1].z += dvH1 * dH1.z;
            dvelocities[j+2].x += dvH2 * dH2.x;
            dvelocities[j+2].y += dvH2 * dH2.y;
            dvelocities[j+2].z += dvH2 * dH2.z;
        }
    }
    
    #ifdef DEBUG 
    printf("q-w: Ecoul = %f Evdw = %f\n", q_energies[0].Ucoul, q_energies[0].Uvdw);
    #endif
}

void calc_nonbonded_qw_forces_host() {
    int mem_size_X = n_atoms_solute * sizeof(coord_t);
    int mem_size_W = 3 * n_waters * sizeof(coord_t);
    int mem_size_DV_X = n_atoms_solute * sizeof(dvel_t);
    int mem_size_DV_W = 3 * n_waters * sizeof(dvel_t);
    int mem_size_MAT = 3 * n_waters * n_qatoms * sizeof(calc_qw_t);

    int n_blocks_q = (n_qatoms + BLOCK_SIZE - 1) / BLOCK_SIZE;
    int n_blocks_w = (n_waters + BLOCK_SIZE - 1) / BLOCK_SIZE;
    //TODO make Evdw & Ecoul work for # of states > 2
    int mem_size_QW_Evdw = min(n_lambdas, 2) * n_blocks_q * n_blocks_w * sizeof(double);
    int mem_size_QW_Ecoul = min(n_lambdas, 2) * n_blocks_q * n_blocks_w * sizeof(double);

    if (A_O == 0) {
        catype_t catype_ow;    // Atom type of first O atom

        catype_ow = catypes[atypes[n_atoms_solute].code - 1];

        A_O = catype_ow.aii_normal;
        B_O = catype_ow.bii_normal;

        #ifdef DEBUG
        printf("Allocating QW_MAT\n");
        #endif
        check_cudaMalloc((void**) &QW_MAT, mem_size_MAT);

        #ifdef DEBUG
        printf("Allocating D_QW_Evdw\n");
        #endif
        check_cudaMalloc((void**) &D_QW_Evdw, mem_size_QW_Evdw);    
        #ifdef DEBUG
        printf("Allocating D_QW_Ecoul\n");
        #endif
        check_cudaMalloc((void**) &D_QW_Ecoul, mem_size_QW_Ecoul);

        check_cudaMalloc((void**) &D_QW_evdw_TOT, sizeof(double));
        check_cudaMalloc((void**) &D_QW_ecoul_TOT, sizeof(double)); 

        h_QW_Evdw = (double*) malloc(mem_size_QW_Evdw);
        h_QW_Ecoul = (double*) malloc(mem_size_QW_Ecoul);

        h_QW_MAT = (calc_qw_t*) malloc(mem_size_MAT);
    }

    cudaMemcpy(X, coords, mem_size_X, cudaMemcpyHostToDevice);
    cudaMemcpy(W, &coords[n_atoms_solute], mem_size_W, cudaMemcpyHostToDevice);
    cudaMemcpy(DV_X, dvelocities, mem_size_DV_X, cudaMemcpyHostToDevice);
    cudaMemcpy(DV_W, &dvelocities[n_atoms_solute], mem_size_DV_W, cudaMemcpyHostToDevice);

    dim3 threads,grid;

    threads = dim3(BLOCK_SIZE, BLOCK_SIZE);
    grid = dim3((n_waters + BLOCK_SIZE - 1) / threads.x, (n_qatoms + BLOCK_SIZE - 1) / threads.y);

    double evdw, ecoul;

    calc_qw_dvel_matrix<<<grid, threads>>>(n_qatoms, n_waters, n_lambdas, crg_ow, crg_hw, A_O, B_O, X, W, D_QW_Evdw, D_QW_Ecoul, QW_MAT, D_qcatypes, D_qatypes, D_qcharges, D_qatoms, D_lambdas, topo);
    calc_qw_dvel_vector_column<<<((n_waters+BLOCK_SIZE - 1) / BLOCK_SIZE), BLOCK_SIZE>>>(n_qatoms, n_waters, DV_X, DV_W, QW_MAT);
    calc_qw_dvel_vector_row<<<((n_qatoms+BLOCK_SIZE - 1) / BLOCK_SIZE), BLOCK_SIZE>>>(n_qatoms, n_waters, DV_X, DV_W, QW_MAT, D_qatoms);

    #ifdef DEBUG
    cudaMemcpy(h_QW_MAT, QW_MAT, mem_size_MAT, cudaMemcpyDeviceToHost);
    #endif

    #ifdef DEBUG
    for (int i = 0; i < n_qatoms; i++) {
        for (int j = 0; j < 3 * n_waters; j++) {
            if (h_QW_MAT[3 * i * n_waters + j].Q.x > 100)
            printf("QW_MAT[%d][%d].Q = %f %f %f\n", i, j, h_QW_MAT[i * 3 * n_waters + j].Q.x, h_QW_MAT[i * 3 * n_waters + j].Q.y, h_QW_MAT[i * 3 * n_waters + j].Q.z);
        }
    }
    #endif

    cudaMemcpy(dvelocities, DV_X, mem_size_DV_X, cudaMemcpyDeviceToHost);
    cudaMemcpy(&dvelocities[n_atoms_solute], DV_W, mem_size_DV_W, cudaMemcpyDeviceToHost);

    //TODO make Evdw & Ecoul work for # of states > 2
    for (int state = 0; state < min(2, n_lambdas); state++) {
        calc_energy_sum<<<1, threads>>>(n_blocks_q, n_blocks_w, D_QW_evdw_TOT, D_QW_ecoul_TOT, &D_QW_Evdw[state * n_blocks_w * n_blocks_q], &D_QW_Ecoul[state * n_blocks_w * n_blocks_q], false);

        cudaMemcpy(&QW_evdw_TOT, D_QW_evdw_TOT, sizeof(double), cudaMemcpyDeviceToHost);
        cudaMemcpy(&QW_ecoul_TOT, D_QW_ecoul_TOT, sizeof(double), cudaMemcpyDeviceToHost);
    
        EQ_nonbond_qw[state].Uvdw += QW_evdw_TOT;
        EQ_nonbond_qw[state].Ucoul += QW_ecoul_TOT;
    }
}

void calc_nonbonded_qq_forces() {
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
        for (int qi = 0; qi < n_qatoms; qi++) {
            for (int qj = qi+1; qj < n_qatoms; qj++) {
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
                    if (q_elscales[k + n_qelscales * state].qi == qi+1 && q_elscales[k + n_qelscales * state].qj == qj+1) {
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
                dva = r2a * ( -Vela -12 * V_a + 6 * V_b) * lambdas[state];

                dvelocities[ai].x -= dva * da.x;
                dvelocities[ai].y -= dva * da.y;
                dvelocities[ai].z -= dva * da.z;

                dvelocities[aj].x += dva * da.x;
                dvelocities[aj].y += dva * da.y;
                dvelocities[aj].z += dva * da.z;
    
                EQ_nonbond_qq[state].Ucoul += Vela;
                EQ_nonbond_qq[state].Uvdw += (V_a - V_b);    
            }
        }
    }

    #ifdef DEBUG 
    printf("q-q: Ecoul = %f Evdw = %f\n", q_energies[0].Ucoul, q_energies[0].Uvdw);
    #endif
}

void calc_qangle_forces(int state) {
    int ic;
    int ai, aj, ak;
    coord_t rji, rjk;
    double bji, bjk;
    double cos_th, th, dth, ener, dv, f1;
    coord_t di, dk;

    for (int i = 0; i < n_qangles; i++) {
        ic = q_angles[i + n_qangles * state].code-1;

        // Skip if angle not present (code 0)
        if (ic == 0) continue;

        ai = q_angles[i + n_qangles * state].ai - 1;
        aj = q_angles[i + n_qangles * state].aj - 1;
        ak = q_angles[i + n_qangles * state].ak - 1;

        rji.x = coords[ai].x - coords[aj].x;
        rji.y = coords[ai].y - coords[aj].y;
        rji.z = coords[ai].z - coords[aj].z;

        rjk.x = coords[ak].x - coords[aj].x;
        rjk.y = coords[ak].y - coords[aj].y;
        rjk.z = coords[ak].z - coords[aj].z;

        bji = sqrt(pow(rji.x, 2) + pow(rji.y, 2) + pow(rji.z, 2));
        bjk = sqrt(pow(rjk.x, 2) + pow(rjk.y, 2) + pow(rjk.z, 2));
        cos_th = rji.x * rjk.x + rji.y * rjk.y + rji.z * rjk.z;
        cos_th /= (bji * bjk);
        if (cos_th > 1) cos_th = 1;
        if (cos_th < -1) cos_th = -1;
        th = acos(cos_th);
        dth = th - to_radians(q_cangles[ic].th0);
        ener = .5 * q_cangles[ic].kth * pow(dth, 2);
        EQ_bond[state].Uangle += ener;

        dv = q_cangles[ic].kth * dth * lambdas[state];
        f1 = sin(th);
        if (abs(f1) < 1E-12) f1 = 1E-12;
        f1 = -1.0 / f1;

        di.x = f1 * (rjk.x / (bji * bjk) - cos_th * rji.x / pow(bji, 2));
        di.y = f1 * (rjk.y / (bji * bjk) - cos_th * rji.y / pow(bji, 2));
        di.z = f1 * (rjk.z / (bji * bjk) - cos_th * rji.z / pow(bji, 2));
        dk.x = f1 * (rji.x / (bji * bjk) - cos_th * rjk.x / pow(bjk, 2));
        dk.y = f1 * (rji.y / (bji * bjk) - cos_th * rjk.y / pow(bjk, 2));
        dk.z = f1 * (rji.z / (bji * bjk) - cos_th * rjk.z / pow(bjk, 2));

        dvelocities[ai].x += dv * di.x;
        dvelocities[ai].y += dv * di.y;
        dvelocities[ai].z += dv * di.z;
        dvelocities[ak].x += dv * dk.x;
        dvelocities[ak].y += dv * dk.y;
        dvelocities[ak].z += dv * dk.z;
        dvelocities[aj].x -= dv * (di.x + dk.x);
        dvelocities[aj].y -= dv * (di.y + dk.y);
        dvelocities[aj].z -= dv * (di.z + dk.z);
    }
}

void calc_qbond_forces(int state) {
    int ic;
    int ai, aj;
    double b, db, ener, dv;
    coord_t rij;

    for (int i = 0; i < n_qbonds; i++) {
        ic = q_bonds[i + n_qbonds * state].code;

        if (ic == 0) continue;

        ai = q_bonds[i + n_qbonds * state].ai - 1;
        aj = q_bonds[i + n_qbonds * state].aj - 1;

        rij.x = coords[aj].x - coords[ai].x;
        rij.y = coords[aj].y - coords[ai].y;
        rij.z = coords[aj].z - coords[ai].z;

        b = sqrt(pow(rij.x, 2) + pow(rij.y, 2) + pow(rij.z, 2));
        db = b - q_cbonds[ic].b0;

        ener = 0.5 * q_cbonds[ic].kb * pow(db, 2);
        EQ_bond[state].Ubond += ener;
        dv = db * q_cbonds[ic].kb * lambdas[state] / b;

        dvelocities[ai].x -= dv * rij.x;
        dvelocities[ai].y -= dv * rij.y;
        dvelocities[ai].z -= dv * rij.z;
        dvelocities[aj].x += dv * rij.x;
        dvelocities[aj].y += dv * rij.y;
        dvelocities[aj].z += dv * rij.z;
    }
}

void calc_qtorsion_forces(int state) {
    int ic;
    int ai, aj, ak, al;
    coord_t rji, rjk, rkl, rnj, rnk, rki, rlj;
    coord_t di, dl, dpi, dpj, dpk, dpl;

    double bj2inv, bk2inv, bjinv, bkinv;
    double bj, bk, cos_phi, phi;
    double arg, dv, f1;
    double ener;

    for (int i = 0; i < n_qtorsions; i++) {
        ic = q_torsions[i + n_qtorsions * state].code;

        if (ic == 0) continue;

        ai = q_torsions[i + n_qtorsions * state].ai - 1;
        aj = q_torsions[i + n_qtorsions * state].aj - 1;
        ak = q_torsions[i + n_qtorsions * state].ak - 1;
        al = q_torsions[i + n_qtorsions * state].al - 1;

        rji.x = coords[ai].x - coords[aj].x;
        rji.y = coords[ai].y - coords[aj].y;
        rji.z = coords[ai].z - coords[aj].z;
        rjk.x = coords[ak].x - coords[aj].x;
        rjk.y = coords[ak].y - coords[aj].y;
        rjk.z = coords[ak].z - coords[aj].z;
        rkl.x = coords[al].x - coords[ak].x;
        rkl.y = coords[al].y - coords[ak].y;
        rkl.z = coords[al].z - coords[ak].z;
        rnj.x = rji.y * rjk.z - rji.z * rjk.y;
        rnj.y = rji.z * rjk.x - rji.x * rjk.z;
        rnj.z = rji.x * rjk.y - rji.y * rjk.x;
        rnk.x = -rjk.y * rkl.z + rjk.z * rkl.y;
        rnk.y = -rjk.z * rkl.x + rjk.x * rkl.z;
        rnk.z = -rjk.x * rkl.y + rjk.y * rkl.x;

        bj = sqrt(pow(rnj.x, 2) + pow(rnj.y, 2) + pow(rnj.z, 2));
        bk = sqrt(pow(rnk.x, 2) + pow(rnk.y, 2) + pow(rnk.z, 2));
        cos_phi = (rnj.x * rnk.x + rnj.y * rnk.y + rnj.z * rnk.z) / (bj * bk);
        if (cos_phi > 1) cos_phi = 1;
        if (cos_phi < -1) cos_phi = -1;
        phi = acos(cos_phi);
        if (rjk.x * (rnj.y * rnk.z - rnj.z * rnk.y)
            + rjk.y * (rnj.z * rnk.x - rnj.x * rnk.z)
            + rjk.z * (rnj.x * rnk.y - rnj.y * rnk.x) < 0) {
            phi = -phi;
        }

        bj2inv = 1 / (pow(rnj.x, 2) + pow(rnj.y, 2) + pow(rnj.z, 2));
        bk2inv = 1 / (pow(rnk.x, 2) + pow(rnk.y, 2) + pow(rnk.z, 2));
        bjinv = sqrt(bj2inv);
        bkinv = sqrt(bk2inv);

        // Energy
        arg = q_ctorsions[ic].n * phi - to_radians(q_ctorsions[ic].d);
        ener = q_ctorsions[ic].k * (1 + cos(arg));
        dv = - q_ctorsions[ic].n * q_ctorsions[ic].k * sin(arg) * lambdas[state];

        // Forces
        f1 = sin(phi);
        if (abs(f1) < 1E-12) f1 = 1E-12;
        f1 = -1 / f1;

        di.x = f1 * (rnk.x * (bjinv * bkinv) - cos_phi * rnj.x * bj2inv);
        di.y = f1 * (rnk.y * (bjinv * bkinv) - cos_phi * rnj.y * bj2inv);
        di.z = f1 * (rnk.z * (bjinv * bkinv) - cos_phi * rnj.z * bj2inv);
        dl.x = f1 * (rnj.x * (bjinv * bkinv) - cos_phi * rnk.x * bk2inv);
        dl.y = f1 * (rnj.y * (bjinv * bkinv) - cos_phi * rnk.y * bk2inv);
        dl.z = f1 * (rnj.z * (bjinv * bkinv) - cos_phi * rnk.z * bk2inv);

        rki.x = rji.x - rjk.x;
        rki.y = rji.y - rjk.y;
        rki.z = rji.z - rjk.z;
        rlj.x = -rjk.x - rkl.x;
        rlj.y = -rjk.y - rkl.y;
        rlj.z = -rjk.z - rkl.z;

        dpi.x = rjk.y * di.z - rjk.z * di.y;
        dpi.y = rjk.z * di.x - rjk.x * di.z;
        dpi.z = rjk.x * di.y - rjk.y * di.x;
        dpj.x = rki.y * di.z - rki.z * di.y + rkl.y * dl.z - rkl.z * dl.y;
        dpj.y = rki.z * di.x - rki.x * di.z + rkl.z * dl.x - rkl.x * dl.z;
        dpj.z = rki.x * di.y - rki.y * di.x + rkl.x * dl.y - rkl.y * dl.x;
        dpk.x = rlj.y * dl.z - rlj.z * dl.y - rji.y * di.z + rji.z * di.y;
        dpk.y = rlj.z * dl.x - rlj.x * dl.z - rji.z * di.x + rji.x * di.z;
        dpk.z = rlj.x * dl.y - rlj.y * dl.x - rji.x * di.y + rji.y * di.x;
        dpl.x = rjk.y * dl.z - rjk.z * dl.y;
        dpl.y = rjk.z * dl.x - rjk.x * dl.z;
        dpl.z = rjk.x * dl.y - rjk.y * dl.x;

        // Update energy and forces
        EQ_bond[state].Utor += ener;

        dvelocities[ai].x += dv * dpi.x;
        dvelocities[ai].y += dv * dpi.y;
        dvelocities[ai].z += dv * dpi.z;

        dvelocities[aj].x += dv * dpj.x;
        dvelocities[aj].y += dv * dpj.y;
        dvelocities[aj].z += dv * dpj.z;

        dvelocities[ak].x += dv * dpk.x;
        dvelocities[ak].y += dv * dpk.y;
        dvelocities[ak].z += dv * dpk.z;

        dvelocities[al].x += dv * dpl.x;
        dvelocities[al].y += dv * dpl.y;
        dvelocities[al].z += dv * dpl.z;
    }
}

/* =============================================
 * == DEVICE
 * =============================================
 */

// Q-W interactions

__device__ void calc_qw_dvel_matrix_incr(int row, int qi, int column, int n_lambdas, int n_qatoms, double crg_ow, double crg_hw, double A_O, double B_O,
    coord_t *Qs, coord_t *Ws, double Evdw_S[BLOCK_SIZE][2 * BLOCK_SIZE], double Ecoul_S[BLOCK_SIZE][2 * BLOCK_SIZE], calc_qw_t *qw,
    q_catype_t *D_qcatypes, q_atype_t *D_qatypes, q_charge_t *D_qcharges, q_atom_t *D_qatoms, double *D_lambdas, topo_t D_topo) {

    int j;
    coord_t dO, dH1, dH2;
    double r2O, rH1, rH2, r6O, rO, r2H1, r2H2;
    double dvO, dvH1, dvH2;
    double V_a, V_b, VelO, VelH1, VelH2;
    q_atype_t qa_type;
    q_catype_t qi_type;
    double ai_aii, ai_bii;

    j = 3 * column;
    dO.x = Ws[j].x - Qs[row].x;
    dO.y = Ws[j].y - Qs[row].y;
    dO.z = Ws[j].z - Qs[row].z;
    dH1.x = Ws[j+1].x - Qs[row].x;
    dH1.y = Ws[j+1].y - Qs[row].y;
    dH1.z = Ws[j+1].z - Qs[row].z;
    dH2.x = Ws[j+2].x - Qs[row].x;
    dH2.y = Ws[j+2].y - Qs[row].y;
    dH2.z = Ws[j+2].z - Qs[row].z;
    
    r2O = pow(dO.x, 2) + pow(dO.y, 2) + pow(dO.z, 2);
    rH1 = sqrt(1.0 / (pow(dH1.x, 2) + pow(dH1.y, 2) + pow(dH1.z, 2)));
    rH2 = sqrt(1.0 / (pow(dH2.x, 2) + pow(dH2.y, 2) + pow(dH2.z, 2)));
    r6O = r2O * r2O * r2O;
    r2O = 1.0 / r2O;
    rO = sqrt(r2O);
    r2H1 = rH1 * rH1;
    r2H2 = rH2 * rH2;

    // Reset potential
    dvO = 0;
    dvH1 = 0;
    dvH2 = 0;

    for (int state = 0; state < n_lambdas; state++) {
        qa_type = D_qatypes[qi + n_qatoms * state];
        qi_type = D_qcatypes[qa_type.code - 1];

        ai_aii = qi_type.Ai;
        ai_bii = qi_type.Bi;

        V_a = ai_aii * A_O / (r6O * r6O);
        V_b = ai_bii * B_O / (r6O);

        VelO = D_topo.coulomb_constant * crg_ow * D_qcharges[qi + n_qatoms * state].q * rO;
        VelH1 = D_topo.coulomb_constant * crg_hw * D_qcharges[qi + n_qatoms * state].q * rH1;
        VelH2 = D_topo.coulomb_constant * crg_hw * D_qcharges[qi + n_qatoms * state].q * rH2;
        
        dvO += r2O * (-VelO - (12 * V_a - 6 * V_b)) * D_lambdas[state];
        dvH1 -= r2H1 * VelH1 * D_lambdas[state];
        dvH2 -= r2H2 * VelH2 * D_lambdas[state];

        // Update Q totals
        Ecoul_S[row][state * BLOCK_SIZE + column] += (VelO + VelH1 + VelH2);
        Evdw_S[row][state * BLOCK_SIZE + column] += (V_a - V_b);
    }

    // Note r6O is not the usual 1/rO^6, but rather rO^6. be careful!!!

    // Update forces on Q-atom
    (*qw).Q.x -= (dvO * dO.x + dvH1 * dH1.x + dvH2 * dH2.x);
    (*qw).Q.y -= (dvO * dO.y + dvH1 * dH1.y + dvH2 * dH2.y);
    (*qw).Q.z -= (dvO * dO.z + dvH1 * dH1.z + dvH2 * dH2.z);

    // Update forces on water
    (*qw).O.x += dvO * dO.x;
    (*qw).O.y += dvO * dO.y;
    (*qw).O.z += dvO * dO.z;
    (*qw).H1.x += dvH1 * dH1.x;
    (*qw).H1.y += dvH1 * dH1.y;
    (*qw).H1.z += dvH1 * dH1.z;
    (*qw).H2.x += dvH2 * dH2.x;
    (*qw).H2.y += dvH2 * dH2.y;
    (*qw).H2.z += dvH2 * dH2.z;
}

__global__ void calc_qw_dvel_matrix(int n_qatoms, int n_waters, int n_lambdas, double crg_ow, double crg_hw, double A_O, double B_O,
    coord_t *X, coord_t *W, double *Evdw, double *Ecoul, calc_qw_t *MAT,
    q_catype_t *D_qcatypes, q_atype_t *D_qatypes, q_charge_t *D_qcharges, q_atom_t *D_qatoms, double *D_lambdas, topo_t D_topo) {
    // Block index
    int bx = blockIdx.x;
    int by = blockIdx.y;

    // Thread index
    int tx = threadIdx.x;
    int ty = threadIdx.y;

    //TODO implement >2 states on GPU
    __shared__ double Evdw_S[BLOCK_SIZE][2 * BLOCK_SIZE];
    __shared__ double Ecoul_S[BLOCK_SIZE][2 * BLOCK_SIZE];    

    Ecoul_S[ty][tx] = 0;
    Evdw_S[ty][tx] = 0;
    Ecoul_S[ty][tx+BLOCK_SIZE] = 0;
    Evdw_S[ty][tx+BLOCK_SIZE] = 0;
    
    int aStart = BLOCK_SIZE * by;
    int bStart = 3 * BLOCK_SIZE * bx;

    if (aStart + ty >= n_qatoms) return;
    if (bStart + 3 * tx >= 3 * n_waters) return;

    __shared__ coord_t Qs[BLOCK_SIZE];
    __shared__ coord_t Ws[3 * BLOCK_SIZE];

    if (tx == 0) {
        Qs[ty] = X[D_qatoms[aStart + ty].a-1];
    }

    if (ty == 0) {
        Ws[3 * tx    ] = W[bStart + 3 * tx    ];
        Ws[3 * tx + 1] = W[bStart + 3 * tx + 1];
        Ws[3 * tx + 2] = W[bStart + 3 * tx + 2];
    }

    __syncthreads();

    calc_qw_t qw;
    memset(&qw, 0, sizeof(calc_qw_t));

    int row = by * BLOCK_SIZE + ty;
    int column = bx * BLOCK_SIZE + tx;
    
    calc_qw_dvel_matrix_incr(ty, aStart + ty, tx, n_lambdas, n_qatoms, crg_ow, crg_hw, A_O, B_O, Qs, Ws, Evdw_S, Ecoul_S,
        &qw, D_qcatypes, D_qatypes, D_qcharges, D_qatoms, D_lambdas, D_topo);

    MAT[column + n_waters * row] = qw;

    __syncthreads();

    int rowlen = (n_waters + BLOCK_SIZE - 1) / BLOCK_SIZE;
    int collen = (n_qatoms + BLOCK_SIZE - 1) / BLOCK_SIZE;

    if (tx == 0 && ty == 0) {
        //TODO implement >2 states on GPU
        for (int state = 0; state < min(2, n_lambdas); state++) {
            double tot_Evdw = 0;
            double tot_Ecoul = 0;
            for (int i = 0; i < BLOCK_SIZE; i++) {
                for (int j = 0; j < BLOCK_SIZE; j++) {
                    tot_Evdw += Evdw_S[i][j + state * BLOCK_SIZE];
                    tot_Ecoul += Ecoul_S[i][j + state * BLOCK_SIZE];
                }
            }
            Evdw[rowlen * collen * state + rowlen * by + bx] = tot_Evdw;
            Ecoul[rowlen * collen * state + rowlen * by + bx] = tot_Ecoul;
        }
    }

    __syncthreads();
}

__global__ void calc_qw_dvel_vector_row(int n_qatoms, int n_waters, dvel_t *DV_X, dvel_t *DV_W, calc_qw_t *MAT, q_atom_t *D_qatoms) {
    int row = blockIdx.x*blockDim.x + threadIdx.x;
    if (row >= n_qatoms) return;

    dvel_t dQ;

    dQ.x = 0;
    dQ.y = 0;
    dQ.z = 0;

    for (int i = 0; i < n_waters; i++) {
        dQ.x += MAT[i + n_waters * row].Q.x;
        dQ.y += MAT[i + n_waters * row].Q.y;
        dQ.z += MAT[i + n_waters * row].Q.z;
    }

    int q = D_qatoms[row].a-1;

    DV_X[q].x += dQ.x;
    DV_X[q].y += dQ.y;
    DV_X[q].z += dQ.z;

    __syncthreads();
}

__global__ void calc_qw_dvel_vector_column(int n_qatoms, int n_waters, dvel_t *DV_X, dvel_t *DV_W, calc_qw_t *MAT) {
    int column = blockIdx.x*blockDim.x + threadIdx.x;
    if (column >= n_waters) return;

    dvel_t dO, dH1, dH2;

    dO.x = 0;
    dO.y = 0;
    dO.z = 0;
    dH1.x = 0;
    dH1.y = 0;
    dH1.z = 0;
    dH2.x = 0;
    dH2.y = 0;
    dH2.z = 0;

    for (int i = 0; i < n_qatoms; i++) {
        dO.x += MAT[column + n_waters * i].O.x;
        dO.y += MAT[column + n_waters * i].O.y;
        dO.z += MAT[column + n_waters * i].O.z;
        dH1.x += MAT[column + n_waters * i].H1.x;
        dH1.y += MAT[column + n_waters * i].H1.y;
        dH1.z += MAT[column + n_waters * i].H1.z;
        dH2.x += MAT[column + n_waters * i].H2.x;
        dH2.y += MAT[column + n_waters * i].H2.y;
        dH2.z += MAT[column + n_waters * i].H2.z;
    }

    DV_W[3*column].x += dO.x;
    DV_W[3*column].y += dO.y;
    DV_W[3*column].z += dO.z;
    DV_W[3*column+1].x += dH1.x;
    DV_W[3*column+1].y += dH1.y;
    DV_W[3*column+1].z += dH1.z;
    DV_W[3*column+2].x += dH2.x;
    DV_W[3*column+2].y += dH2.y;
    DV_W[3*column+2].z += dH2.z;

    __syncthreads();
}

// Q-Q interactions

// Q-P interactions

__device__ void calc_qp_dvel_matrix_incr(int row, int qi, int column, int pj, int n_lambdas, int n_qatoms,
    coord_t *Qs, coord_t *Ps, int *LJs, bool *excluded_s, double Evdw_S[BLOCK_SIZE][2 * BLOCK_SIZE], double Ecoul_S[BLOCK_SIZE][2 * BLOCK_SIZE], calc_qp_t *qp,
    q_catype_t *D_qcatypes, q_atype_t *D_qatypes, q_charge_t *D_qcharges, p_atom_t *D_patoms, q_atom_t *D_qatoms, double *D_lambdas,
    catype_t *D_catypes, atype_t *D_atypes, ccharge_t *D_ccharges, charge_t *D_charges, topo_t D_topo) {
    
    coord_t da;
    double r2, r6, r;
    double ai_aii, aj_aii, ai_bii, aj_bii;
    q_catype_t qi_type;
    catype_t aj_type;
    bool bond23, bond14;
    double scaling, Vel, V_a, V_b, dv;

    int j = D_patoms[pj].a-1;

    bond23 = LJs[row * BLOCK_SIZE + column] == 3;
    bond14 = LJs[row * BLOCK_SIZE + column] == 1;

    if (bond23) return;
    if (excluded_s[row] || excluded_s[BLOCK_SIZE + column]) return;

    scaling = bond14 ? D_topo.el14_scale : 1;

    da.x = Qs[row].x - Ps[column].x;
    da.y = Qs[row].y - Ps[column].y;
    da.z = Qs[row].z - Ps[column].z;

    r2 = pow(da.x, 2) + pow(da.y, 2) + pow(da.z, 2);

    r6 = r2 * r2 * r2;
    r2 = 1 / r2;
    r = sqrt(r2);

    for (int state = 0; state < n_lambdas; state++) {
        qi_type = D_qcatypes[D_qatypes[qi + n_qatoms * state].code - 1];
        aj_type = D_catypes[D_atypes[j].code - 1];

        ai_aii = bond14 ? qi_type.Ai_14 : qi_type.Ai;
        aj_aii = bond14 ? aj_type.aii_1_4 : aj_type.aii_normal;
        ai_bii = bond14 ? qi_type.Bi_14 : qi_type.Bi;
        aj_bii = bond14 ? aj_type.bii_1_4 : aj_type.bii_normal;

        Vel = D_topo.coulomb_constant * scaling * D_qcharges[qi + n_qatoms * state].q * D_ccharges[D_charges[j].code - 1].charge * r;
        V_a = ai_aii * aj_aii / (r6 * r6);
        V_b = ai_bii * aj_bii / r6;
        dv = r2 * (-Vel - (12 * V_a - 6 * V_b)) * D_lambdas[state];

        // if (state == 0 && qi == 0 && pj == 1) {
        //     printf("crg_q = %f crg_j = %f r = %f\n", D_qcharges[qi + n_qatoms * state].q, D_ccharges[D_charges[pj].code - 1].charge, r);
        //     printf("ai_aii = %f aj_aii = %f ai_bii = %f aj_bii = %f\n", ai_aii, aj_aii, ai_bii, aj_bii);
        // }

        // Update forces
        qp->Q.x += dv * da.x;
        qp->Q.y += dv * da.y;
        qp->Q.z += dv * da.z;
        qp->P.x -= dv * da.x;
        qp->P.y -= dv * da.y;
        qp->P.z -= dv * da.z;

        // Update Q totals
        Ecoul_S[row][state * BLOCK_SIZE + column] += Vel;
        Evdw_S[row][state * BLOCK_SIZE + column] += (V_a - V_b);
    }
}

__global__ void calc_qp_dvel_matrix(int n_qatoms, int n_patoms, int n_lambdas, int n_atoms_solute,
    coord_t *X, double *Evdw, double *Ecoul, calc_qp_t *QP_MAT,
    q_catype_t *D_qcatypes, q_atype_t *D_qatypes, q_charge_t *D_qcharges, p_atom_t *D_patoms, q_atom_t *D_qatoms, double *D_lambdas,
    int *D_LJ_matrix, bool *D_excluded, catype_t *D_catypes, atype_t *D_atypes, ccharge_t *D_ccharges, charge_t *D_charges, topo_t D_topo) {
    
    // Block index
    int bx = blockIdx.x;
    int by = blockIdx.y;

    // Thread index
    int tx = threadIdx.x;
    int ty = threadIdx.y;

    //TODO implement >2 states on GPU
    __shared__ double Evdw_S[BLOCK_SIZE][2 * BLOCK_SIZE];
    __shared__ double Ecoul_S[BLOCK_SIZE][2 * BLOCK_SIZE];    

    Ecoul_S[ty][tx] = 0;
    Evdw_S[ty][tx] = 0;
    Ecoul_S[ty][tx+BLOCK_SIZE] = 0;
    Evdw_S[ty][tx+BLOCK_SIZE] = 0;

    int aStart = BLOCK_SIZE * by;
    int bStart = BLOCK_SIZE * bx;

    if (aStart + ty >= n_qatoms) return;
    if (bStart + tx >= n_patoms) return;

    int qi = D_qatoms[aStart + ty].a-1;
    int pj = D_patoms[bStart + tx].a-1;

    __shared__ coord_t Qs[BLOCK_SIZE];
    __shared__ coord_t Ps[BLOCK_SIZE];
    __shared__     int LJs[BLOCK_SIZE * BLOCK_SIZE];
    __shared__    bool excluded_s[2 * BLOCK_SIZE];

    if (tx == 0) {
        Qs[ty] = X[qi];
        excluded_s[ty] = D_excluded[qi];
    }

    if (ty == 0) {
        Ps[tx] = X[pj];
        excluded_s[BLOCK_SIZE + tx] = D_excluded[pj];
    }
    LJs[ty * BLOCK_SIZE + tx] = D_LJ_matrix[qi * n_atoms_solute + pj];

    __syncthreads();

    calc_qp_t qp;
    memset(&qp, 0, sizeof(calc_qw_t));

    int row = by * BLOCK_SIZE + ty;
    int column = bx * BLOCK_SIZE + tx;
    
    calc_qp_dvel_matrix_incr(ty, row, tx, column, n_lambdas, n_qatoms, Qs, Ps, LJs, excluded_s, Evdw_S, Ecoul_S,
        &qp, D_qcatypes, D_qatypes, D_qcharges, D_patoms, D_qatoms, D_lambdas, D_catypes, D_atypes, D_ccharges, D_charges, D_topo);

    QP_MAT[n_patoms * row + column] = qp;

    __syncthreads();

    int rowlen = (n_patoms + BLOCK_SIZE - 1) / BLOCK_SIZE;
    int collen = (n_qatoms + BLOCK_SIZE - 1) / BLOCK_SIZE;

    if (tx == 0 && ty == 0) {
        //TODO implement >2 states on GPU
        for (int state = 0; state < min(2, n_lambdas); state++) {
            double tot_Evdw = 0;
            double tot_Ecoul = 0;
            for (int i = 0; i < BLOCK_SIZE; i++) {
                for (int j = 0; j < BLOCK_SIZE; j++) {
                    tot_Evdw += Evdw_S[i][j + state * BLOCK_SIZE];
                    tot_Ecoul += Ecoul_S[i][j + state * BLOCK_SIZE];
                }
            }
            Evdw[rowlen * collen * state + rowlen * by + bx] = tot_Evdw;
            Ecoul[rowlen * collen * state + rowlen * by + bx] = tot_Ecoul;
        }
    }

    __syncthreads();
}

__global__ void calc_qp_dvel_vector_row(int n_qatoms, int n_patoms, dvel_t *DV_X, calc_qp_t *QP_MAT, q_atom_t *D_qatoms) {
    int row = blockIdx.x*blockDim.x + threadIdx.x;
    if (row >= n_qatoms) return;

    dvel_t dQ;

    dQ.x = 0;
    dQ.y = 0;
    dQ.z = 0;

    for (int i = 0; i < n_patoms; i++) {
        dQ.x += QP_MAT[i + n_patoms * row].Q.x;
        dQ.y += QP_MAT[i + n_patoms * row].Q.y;
        dQ.z += QP_MAT[i + n_patoms * row].Q.z;
    }

    int q = D_qatoms[row].a-1;

    DV_X[q].x += dQ.x;
    DV_X[q].y += dQ.y;
    DV_X[q].z += dQ.z;

    __syncthreads();
}

__global__ void calc_qp_dvel_vector_column(int n_qatoms, int n_patoms, dvel_t *DV_X, calc_qp_t *QP_MAT, p_atom_t *D_patoms) {
    int column = blockIdx.x*blockDim.x + threadIdx.x;
    if (column >= n_patoms) return;

    dvel_t dP;

    dP.x = 0;
    dP.y = 0;
    dP.z = 0;

    for (int i = 0; i < n_qatoms; i++) {
        dP.x += QP_MAT[column + n_patoms * i].P.x;
        dP.y += QP_MAT[column + n_patoms * i].P.y;
        dP.z += QP_MAT[column + n_patoms * i].P.z;
    }

    int p = D_patoms[column].a-1;

    DV_X[p].x += dP.x;
    DV_X[p].y += dP.y;
    DV_X[p].z += dP.z;

    __syncthreads();
}

void clean_d_qatoms() {
    cudaFree(QW_MAT);
    cudaFree(QP_MAT);
    cudaFree(QQ_MAT);
    cudaFree(D_qcatypes);
    cudaFree(D_qatypes);
    cudaFree(D_qcharges);
    cudaFree(D_qatoms);
    cudaFree(D_lambdas);
    free(h_QW_MAT);
    free(h_QP_MAT);
}
