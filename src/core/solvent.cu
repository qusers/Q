#include "system.h"
#include "solvent.h"
#include "utils.h"

#include <stdio.h>
#include <time.h>
#include <unistd.h>

/* =============================================
 * == SOLVENT INTERACTIONS
 * =============================================
 */

coord_t *P, *X;
dvel_t *DV, *DV_P;
calc_ww_t *WW_MAT, *h_WW_MAT;
calc_pw_t *PW_MAT, *h_PW_MAT;

// Constants pointers
ccharge_t *D_ccharges;
charge_t *D_charges;
catype_t *D_catypes;
atype_t * D_atypes;
p_atom_t *D_patoms;
bool gpu_set = false;

//ONLY call if there are actually solvent atoms, or get segfaulted
void calc_nonbonded_ww_forces() {
    double rOX, rH1X, rH2X, r2;
    coord_t dOX, dH1X, dH2X;
    double Vel, V_a, V_b, dv;

    // Initialize water constants
    if (A_OO == 0) {
        catype_t catype_ow;    // Atom type of first O, H atom
        ccharge_t ccharge_ow, ccharge_hw; // Charge of first O, H atom

        catype_ow = catypes[atypes[n_atoms_solute].code - 1];
        ccharge_ow = ccharges[charges[n_atoms_solute].code - 1];
        ccharge_hw = ccharges[charges[n_atoms_solute+1].code - 1];

        A_OO = pow(catype_ow.aii_normal, 2);
        B_OO = pow(catype_ow.bii_normal, 2);

        crg_ow = ccharge_ow.charge;
        crg_hw = ccharge_hw.charge;
    }

    for (int i = n_atoms_solute; i < n_atoms; i+=3) {
        for (int j = i+3; j < n_atoms; j+=3) {
            if (excluded[i] || excluded[j]) continue;
            // --- O - (O,H1,H2) ---
            dOX.x = coords[j].x - coords[i].x;
            dOX.y = coords[j].y - coords[i].y;
            dOX.z = coords[j].z - coords[i].z;
            rOX = pow(dOX.x, 2) + pow(dOX.y, 2) + pow(dOX.z, 2);

            // O - H1
            j += 1;

            dH1X.x = coords[j].x - coords[i].x;
            dH1X.y = coords[j].y - coords[i].y;
            dH1X.z = coords[j].z - coords[i].z;
            rH1X = pow(dH1X.x, 2) + pow(dH1X.y, 2) + pow(dH1X.z, 2);

            // O-H2 (X=H2)
            j += 1;

            dH2X.x = coords[j].x - coords[i].x;
            dH2X.y = coords[j].y - coords[i].y;
            dH2X.z = coords[j].z - coords[i].z;
            rH2X = pow(dH2X.x, 2) + pow(dH2X.y, 2) + pow(dH2X.z, 2);
            rOX = sqrt(1 / rOX);
            rH1X = sqrt(1 / rH1X);
            rH2X = sqrt(1 / rH2X);

            // O - O
            r2 = rOX * rOX;
            Vel = Coul * pow(crg_ow, 2) * rOX;
            V_a = A_OO * (r2*r2*r2) * (r2*r2*r2);
            V_b = B_OO * (r2*r2*r2);
            energies.Uvdw += (V_a - V_b);
            energies.Ucoul += Vel;
            dv = r2 * (-Vel - 12 * V_a + 6 * V_b);
            j -= 2; //move pointer back to O in interacting molecule
            dvelocities[i].x -= (dv * dOX.x);
            dvelocities[j].x += (dv * dOX.x);
            dvelocities[i].y -= (dv * dOX.y);
            dvelocities[j].y += (dv * dOX.y);
            dvelocities[i].z -= (dv * dOX.z);
            dvelocities[j].z += (dv * dOX.z);

            // O - H1
            r2 = pow(rH1X, 2);
            Vel = Coul * crg_ow * crg_hw * rH1X;
            energies.Ucoul += Vel;
            dv = r2 * (-Vel);
            j += 1; //point to H1 in j-molecule
            dvelocities[i].x -= (dv * dH1X.x);
            dvelocities[j].x += (dv * dH1X.x);
            dvelocities[i].y -= (dv * dH1X.y);
            dvelocities[j].y += (dv * dH1X.y);
            dvelocities[i].z -= (dv * dH1X.z);
            dvelocities[j].z += (dv * dH1X.z);

            // O - H2
            r2 = pow(rH2X, 2);
            Vel = Coul * crg_ow * crg_hw * rH2X;
            energies.Ucoul += Vel;
            dv = r2 * (-Vel);
            j += 1; //point to H2 in j-molecule
            dvelocities[i].x -= (dv * dH2X.x);
            dvelocities[j].x += (dv * dH2X.x);
            dvelocities[i].y -= (dv * dH2X.y);
            dvelocities[j].y += (dv * dH2X.y);
            dvelocities[i].z -= (dv * dH2X.z);
            dvelocities[j].z += (dv * dH2X.z);

            // --- H1 - (O,H1,H2) ---
            i += 1; //Point to H1 in i-molecule
            j -= 2; //Point to O in j-molecule

            // H1 - O (X=O)
            dOX.x = coords[j].x - coords[i].x;
            dOX.y = coords[j].y - coords[i].y;
            dOX.z = coords[j].z - coords[i].z;
            rOX = pow(dOX.x, 2) + pow(dOX.y, 2) + pow(dOX.z, 2);

            // H1 - H1 (X=H1)
            j += 1;

            dH1X.x = coords[j].x - coords[i].x;
            dH1X.y = coords[j].y - coords[i].y;
            dH1X.z = coords[j].z - coords[i].z;
            rH1X = pow(dH1X.x, 2) + pow(dH1X.y, 2) + pow(dH1X.z, 2);

            // H1 - H2 (X=H2)
            j += 1;

            dH2X.x = coords[j].x - coords[i].x;
            dH2X.y = coords[j].y - coords[i].y;
            dH2X.z = coords[j].z - coords[i].z;
            rH2X = pow(dH2X.x, 2) + pow(dH2X.y, 2) + pow(dH2X.z, 2);
            rOX = sqrt(1 / rOX);
            rH1X = sqrt(1 / rH1X);
            rH2X = sqrt(1 / rH2X);

            // H1 - O
            r2 = rOX * rOX;
            Vel = Coul * crg_hw * crg_ow * rOX;
            energies.Ucoul += Vel;
            dv = r2 * (-Vel);
            j -= 2; //move pointer back to O in interacting molecule
            dvelocities[i].x -= (dv * dOX.x);
            dvelocities[j].x += (dv * dOX.x);
            dvelocities[i].y -= (dv * dOX.y);
            dvelocities[j].y += (dv * dOX.y);
            dvelocities[i].z -= (dv * dOX.z);
            dvelocities[j].z += (dv * dOX.z);

            // H1 - H1
            r2 = pow(rH1X, 2);
            Vel = Coul * crg_hw * crg_hw * rH1X;
            energies.Ucoul += Vel;
            dv = r2 * (-Vel);
            j += 1; //point to H1 in j-molecule
            dvelocities[i].x -= (dv * dH1X.x);
            dvelocities[j].x += (dv * dH1X.x);
            dvelocities[i].y -= (dv * dH1X.y);
            dvelocities[j].y += (dv * dH1X.y);
            dvelocities[i].z -= (dv * dH1X.z);
            dvelocities[j].z += (dv * dH1X.z);

            // H1 - H2
            r2 = pow(rH2X, 2);
            Vel = Coul * crg_hw * crg_hw * rH2X;
            energies.Ucoul += Vel;
            dv = r2 * (-Vel);
            j += 1; //point to H2 in j-molecule
            dvelocities[i].x -= (dv * dH2X.x);
            dvelocities[j].x += (dv * dH2X.x);
            dvelocities[i].y -= (dv * dH2X.y);
            dvelocities[j].y += (dv * dH2X.y);
            dvelocities[i].z -= (dv * dH2X.z);
            dvelocities[j].z += (dv * dH2X.z);

            // --- H2 - (O,H1,H2) ---
            i += 1; //Point to H2 in i-molecule
            j -= 2; //Point to O in j-molecule

            // H2 - O (X=O)
            dOX.x = coords[j].x - coords[i].x;
            dOX.y = coords[j].y - coords[i].y;
            dOX.z = coords[j].z - coords[i].z;
            rOX = pow(dOX.x, 2) + pow(dOX.y, 2) + pow(dOX.z, 2);

            // H2 - H1 (X=H1)
            j += 1;

            dH1X.x = coords[j].x - coords[i].x;
            dH1X.y = coords[j].y - coords[i].y;
            dH1X.z = coords[j].z - coords[i].z;
            rH1X = pow(dH1X.x, 2) + pow(dH1X.y, 2) + pow(dH1X.z, 2);

            // H2 - H2 (X=H2)
            j += 1;

            dH2X.x = coords[j].x - coords[i].x;
            dH2X.y = coords[j].y - coords[i].y;
            dH2X.z = coords[j].z - coords[i].z;
            rH2X = pow(dH2X.x, 2) + pow(dH2X.y, 2) + pow(dH2X.z, 2);
            rOX = sqrt(1 / rOX);
            rH1X = sqrt(1 / rH1X);
            rH2X = sqrt(1 / rH2X);

            // H2 - O
            r2 = rOX * rOX;
            Vel = Coul * crg_hw * crg_ow * rOX;
            energies.Ucoul += Vel;
            dv = r2 * (-Vel);
            j -= 2; //move pointer back to O in interacting molecule
            dvelocities[i].x -= (dv * dOX.x);
            dvelocities[j].x += (dv * dOX.x);
            dvelocities[i].y -= (dv * dOX.y);
            dvelocities[j].y += (dv * dOX.y);
            dvelocities[i].z -= (dv * dOX.z);
            dvelocities[j].z += (dv * dOX.z);

            // H2 - H1
            r2 = pow(rH1X, 2);
            Vel = Coul * crg_hw * crg_hw * rH1X;
            energies.Ucoul += Vel;
            dv = r2 * (-Vel);
            j += 1; //point to H1 in j-molecule
            dvelocities[i].x -= (dv * dH1X.x);
            dvelocities[j].x += (dv * dH1X.x);
            dvelocities[i].y -= (dv * dH1X.y);
            dvelocities[j].y += (dv * dH1X.y);
            dvelocities[i].z -= (dv * dH1X.z);
            dvelocities[j].z += (dv * dH1X.z);

            // H1 - H2
            r2 = pow(rH2X, 2);
            Vel = Coul * crg_hw * crg_hw * rH2X;
            energies.Ucoul += Vel;
            dv = r2 * (-Vel);
            j += 1; //point to H2 in j-molecule
            dvelocities[i].x -= (dv * dH2X.x);
            dvelocities[j].x += (dv * dH2X.x);
            dvelocities[i].y -= (dv * dH2X.y);
            dvelocities[j].y += (dv * dH2X.y);
            dvelocities[i].z -= (dv * dH2X.z);
            dvelocities[j].z += (dv * dH2X.z);

            // Put i and j indexes back to prevent troubles
            i -= 2;
            j -= 2;
        }
    }

    // for (int i = 0; i < n_atoms; i++) {
    //     printf("dvelocities[%d] = %f %f %f\n", i, dvelocities[i].x, dvelocities[i].y, dvelocities[i].z);
    // }

    #ifdef DEBUG 
    printf("solvent: Ecoul = %f Evdw = %f\n", energies.Ucoul, energies.Uvdw);
    #endif
}

void calc_nonbonded_ww_forces_host() {
    int mem_size_X = 3 * n_waters * sizeof(coord_t);
    int mem_size_MAT = n_waters * n_waters * sizeof(calc_ww_t);
    int mem_size_DV = 3 * n_waters * sizeof(dvel_t);

    cudaError_t error;

    // Initialize water constants
    if (A_OO == 0) {
        catype_t catype_ow;    // Atom type of first O, H atom
        ccharge_t ccharge_ow, ccharge_hw; // Charge of first O, H atom

        catype_ow = catypes[atypes[n_atoms_solute].code - 1];
        ccharge_ow = ccharges[charges[n_atoms_solute].code - 1];
        ccharge_hw = ccharges[charges[n_atoms_solute+1].code - 1];

        A_OO = pow(catype_ow.aii_normal, 2);
        B_OO = pow(catype_ow.bii_normal, 2);

        crg_ow = ccharge_ow.charge;
        crg_hw = ccharge_hw.charge;

        #ifdef DEBUG
        printf("Allocating X\n");
        #endif
        check_cudaMalloc((void**) &X, mem_size_X);
        #ifdef DEBUG
        printf("Allocating WW_MAT\n");
        #endif
        check_cudaMalloc((void**) &WW_MAT, mem_size_MAT);
        #ifdef DEBUG
        printf("Allocating DV\n");
        #endif
        check_cudaMalloc((void**) &DV, mem_size_DV);    
        h_WW_MAT = (calc_ww_t*) malloc(mem_size_MAT);
    }

    cudaMemcpy(X, &coords[n_atoms_solute], mem_size_X, cudaMemcpyHostToDevice);
    cudaMemset(WW_MAT, 0, mem_size_MAT);
    cudaMemcpy(DV, &dvelocities[n_atoms_solute], mem_size_DV, cudaMemcpyHostToDevice);

    dim3 threads,grid;

    threads = dim3(BLOCK_SIZE, BLOCK_SIZE);
    grid = dim3((n_waters + BLOCK_SIZE - 1) / threads.x, (n_waters + BLOCK_SIZE - 1) / threads.y);

    double evdw, ecoul;

    calc_ww_dvel_matrix<<<grid, threads>>>(n_waters, crg_ow, crg_hw, A_OO, B_OO, X, &evdw, &ecoul, WW_MAT);
    calc_ww_dvel_vector<<<((n_waters+BLOCK_SIZE - 1) / BLOCK_SIZE), BLOCK_SIZE>>>(n_waters, DV, WW_MAT);

    // cudaMemcpy(h_WW_MAT, WW_MAT, mem_size_MAT, cudaMemcpyDeviceToHost);

    // for (int i = 0; i < n_waters; i++) {
    //     printf("X[%d] = %f %f %f\n", i, coords[i].x, coords[i].y, coords[i].z);
    // }

    // for (int i = 0; i < n_waters; i++) {
    //     for (int j = 0; j < n_waters; j++) {
    //         printf("MAT[%d][%d].O = %f %f %f\n", i, j, h_MAT[i * n_waters + j].O.x, h_MAT[i * n_waters + j].O.y, h_MAT[i * n_waters + j].O.z);
    //     }
    // }

    cudaMemcpy(&dvelocities[n_atoms_solute], DV, mem_size_DV, cudaMemcpyDeviceToHost);

    // for (int i = 0; i < n_waters; i++) {
    //     printf("dvelocities[%d] = %f %f %f\n", i, dvelocities[i].x, dvelocities[i].y, dvelocities[i].z);
    // }
}


void calc_nonbonded_pw_forces() {
    coord_t da;
    double r2a, ra, r6a;
    double Vela, V_a, V_b;
    double dva;
    double qi, qj;
    double ai_aii, aj_aii, ai_bii, aj_bii;
    catype_t ai_type, aj_type;
    int i;

    for (int pi = 0; pi < n_patoms; pi++) {
        for (int j = n_atoms_solute; j < n_atoms; j++) {
            i = p_atoms[pi].a-1;
            // if (excluded[i] || excluded[j]) continue;
            qi = ccharges[charges[i].code - 1].charge;
            qj = ccharges[charges[j].code - 1].charge;

            ai_type = catypes[atypes[i].code - 1];
            aj_type = catypes[atypes[j].code - 1];

            da.x = coords[j].x - coords[i].x;
            da.y = coords[j].y - coords[i].y;
            da.z = coords[j].z - coords[i].z;
            r2a = 1 / (pow(da.x, 2) + pow(da.y, 2) + pow(da.z, 2));
            ra = sqrt(r2a);
            r6a = r2a * r2a * r2a;

            Vela = Coul * qi * qj * ra;

            ai_aii = ai_type.aii_normal;
            aj_aii = aj_type.aii_normal;
            ai_bii = ai_type.bii_normal;
            aj_bii = aj_type.bii_normal;

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

            // if (j - n_atoms_solute == 1754) {
            //     printf("MAT[%d][%d].W = %f %f %f\n", i, j - n_atoms_solute, dva * da.x, dva * da.y, dva * da.z);
            // }
        }
    }

    // for (int i = 0; i < n_atoms; i++) {
    //     printf("dvelocities[%d] = %f %f %f\n", i, dvelocities[i].x, dvelocities[i].y, dvelocities[i].z);
    // }

    #ifdef DEBUG 
    printf("solute-solvent: Ecoul = %f Evdw = %f\n", energies.Ucoul, energies.Uvdw);
    #endif 
}

void calc_nonbonded_pw_forces_host() {
    int mem_size_X = 3 * n_waters * sizeof(coord_t);
    int mem_size_P = n_patoms * sizeof(coord_t);
    int mem_size_DV = 3 * n_waters * sizeof(dvel_t);
    int mem_size_DV_P = n_patoms * sizeof(dvel_t);
    int mem_size_PW_MAT = 3 * n_waters * n_patoms * sizeof(calc_pw_t);

    int mem_size_ccharges = n_ccharges * sizeof(ccharge_t);
    int mem_size_charges = n_atoms * sizeof(charge_t);
    int mem_size_catypes = n_catypes * sizeof(catype_t);
    int mem_size_atypes = n_atoms * sizeof(atype_t);
    int mem_size_patoms = n_patoms * sizeof(p_atom_t);

    if (!gpu_set) {
        gpu_set = true;

        #ifdef DEBUG
        printf("Allocating P\n");
        #endif
        check_cudaMalloc((void**) &P, mem_size_P);
        #ifdef DEBUG
        printf("Allocating DV_P\n");
        #endif
        check_cudaMalloc((void**) &DV_P, mem_size_DV_P);
        #ifdef DEBUG
        printf("Allocating PW_MAT\n");
        #endif
        check_cudaMalloc((void**) &PW_MAT, mem_size_PW_MAT);

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
        printf("All GPU solvent memory allocated\n");
        #endif

        cudaMemcpy(D_ccharges, ccharges, mem_size_ccharges, cudaMemcpyHostToDevice);
        cudaMemcpy(D_charges, charges, mem_size_charges, cudaMemcpyHostToDevice);
        cudaMemcpy(D_catypes, catypes, mem_size_catypes, cudaMemcpyHostToDevice);
        cudaMemcpy(D_atypes, atypes, mem_size_atypes, cudaMemcpyHostToDevice);
        cudaMemcpy(D_patoms, p_atoms, mem_size_patoms, cudaMemcpyHostToDevice);

        h_PW_MAT = (calc_pw_t*) malloc(mem_size_PW_MAT);
    }

    cudaMemcpy(P, coords, mem_size_P, cudaMemcpyHostToDevice);
    cudaMemcpy(X, &coords[n_atoms_solute], mem_size_X, cudaMemcpyHostToDevice);
    cudaMemcpy(DV_P, dvelocities, mem_size_DV_P, cudaMemcpyHostToDevice);
    cudaMemcpy(DV, &dvelocities[n_atoms_solute], mem_size_DV, cudaMemcpyHostToDevice);

    dim3 threads,grid;

    threads = dim3(BLOCK_SIZE, BLOCK_SIZE);
    grid = dim3((3*n_waters + BLOCK_SIZE - 1) / threads.x, (n_patoms + BLOCK_SIZE - 1) / threads.y);
    double evdw, ecoul;

    // __global__ void calc_pw_dvel_matrix(int n_patoms, int n_waters,
    //     coord_t *P, coord_t *X, double *Evdw, double *Ecoul, calc_pw_t *PW_MAT,
    //     ccharge_t *D_ccharges, charge_t *D_charges, catype_t *D_catypes, atype_t *D_atypes, p_atom_t *D_patoms);
    // __global__ void calc_pw_dvel_vector_row(int n_patoms, int n_waters, dvel_t *DV_P, dvel_t *DV, calc_pw_t *PW_MAT);
    // __global__ void calc_pw_dvel_vector_column(int n_patoms, int n_waters, dvel_t *DV_P, dvel_t *DV, calc_pw_t *PW_MAT);
    
    calc_pw_dvel_matrix<<<grid, threads>>>(n_patoms, n_waters, P, X, &evdw, &ecoul, PW_MAT, D_ccharges, D_charges, D_catypes, D_atypes, D_patoms);
    calc_pw_dvel_vector_column<<<((3*n_waters+BLOCK_SIZE - 1) / BLOCK_SIZE), BLOCK_SIZE>>>(n_patoms, n_waters, DV_P, DV, PW_MAT);
    calc_pw_dvel_vector_row<<<((n_patoms+BLOCK_SIZE - 1) / BLOCK_SIZE), BLOCK_SIZE>>>(n_patoms, n_waters, DV_P, DV, PW_MAT);

    cudaMemcpy(h_PW_MAT, PW_MAT, mem_size_PW_MAT, cudaMemcpyDeviceToHost);

    // for (int i = 0; i < n_waters; i++) {
    //     printf("X[%d] = %f %f %f\n", i, coords[i].x, coords[i].y, coords[i].z);
    // }

    // printf("n_patoms = %d n_watoms = %d\n", n_patoms, 3 * n_waters);
    // for (int i = 0; i < n_patoms; i++) {
        // for (int j = 1745; j < 3 * n_waters; j++) {
            // printf("MAT[%d][%d].W = %f %f %f\n", i, 1754, h_PW_MAT[3 * i * n_waters + 1754].W.x, h_PW_MAT[3 * i * n_waters + 1754].W.y, h_PW_MAT[3 * i * n_waters + 1754].W.z);
        // }
    // }

    cudaMemcpy(dvelocities, DV_P, mem_size_DV_P, cudaMemcpyDeviceToHost);
    cudaMemcpy(&dvelocities[n_atoms_solute], DV, mem_size_DV, cudaMemcpyDeviceToHost);

    // for (int i = 0; i < n_atoms; i++) {
    //     printf("dvelocities[%d] = %f %f %f\n", i, dvelocities[i].x, dvelocities[i].y, dvelocities[i].z);
    // }
}

/* =============================================
 * == DEVICE
 * =============================================
 */

__device__ void set_water(int n_waters, int row, int column, dvel_t *val, calc_ww_t *MAT) {
    // if (row < 32 && column < 32) {
    //     printf("MAT[%d][%d].O = %f %f %f\n", row, column, val[0].x, val[0].y, val[0].z);
    // }
    MAT[column + n_waters * row].O  = val[0];
    MAT[column + n_waters * row].H1 = val[1];
    MAT[column + n_waters * row].H2 = val[2];
}

__device__ void calc_ww_dvel_matrix_incr(int row, int column, double crg_ow, double crg_hw, double A_OO, double B_OO,
    coord_t *Xs, coord_t *Ys, double *Evdw, double *Ecoul, dvel_t *water_a, dvel_t *water_b) {
    double rOX, rH1X, rH2X, r2;
    coord_t dOX, dH1X, dH2X;
    double Vel, V_a, V_b, dv;
    double tempX, tempY, tempZ;
    
    int i = 3 * row;
    int j = 3 * column;
    int wi = 0, wj = 0;

    // --- O - (O,H1,H2) ---
    dOX.x = Ys[j].x - Xs[i].x;
    dOX.y = Ys[j].y - Xs[i].y;
    dOX.z = Ys[j].z - Xs[i].z;
    rOX = pow(dOX.x, 2) + pow(dOX.y, 2) + pow(dOX.z, 2);

    // O - H1
    j += 1;
    wj += 1;

    dH1X.x = Ys[j].x - Xs[i].x;
    dH1X.y = Ys[j].y - Xs[i].y;
    dH1X.z = Ys[j].z - Xs[i].z;
    rH1X = pow(dH1X.x, 2) + pow(dH1X.y, 2) + pow(dH1X.z, 2);

    // O-H2 (X=H2)
    j += 1;
    wj += 1;

    dH2X.x = Ys[j].x - Xs[i].x;
    dH2X.y = Ys[j].y - Xs[i].y;
    dH2X.z = Ys[j].z - Xs[i].z;
    rH2X = pow(dH2X.x, 2) + pow(dH2X.y, 2) + pow(dH2X.z, 2);
    rOX = sqrt(1 / rOX);
    rH1X = sqrt(1 / rH1X);
    rH2X = sqrt(1 / rH2X);

    // O - O
    r2 = rOX * rOX;
    Vel = Coul * pow(crg_ow, 2) * rOX;
    V_a = A_OO * (r2*r2*r2) * (r2*r2*r2);
    V_b = B_OO * (r2*r2*r2);
    *Evdw += (V_a - V_b);
    *Ecoul += Vel;
    dv = r2 * (-Vel - 12 * V_a + 6 * V_b);
    j -= 2; //move pointer back to O in interacting molecule
    wj -= 2;

    water_a[wi].x -= (dv * dOX.x);
    water_b[wj].x += (dv * dOX.x);
    water_a[wi].y -= (dv * dOX.y);
    water_b[wj].y += (dv * dOX.y);
    water_a[wi].z -= (dv * dOX.z);
    water_b[wj].z += (dv * dOX.z);

    // O - H1
    r2 = pow(rH1X, 2);
    Vel = Coul * crg_ow * crg_hw * rH1X;
    *Ecoul += Vel;
    dv = r2 * (-Vel);
    j += 1; //point to H1 in j-molecule
    wj += 1;

    water_a[wi].x -= (dv * dH1X.x);
    water_b[wj].x += (dv * dH1X.x);
    water_a[wi].y -= (dv * dH1X.y);
    water_b[wj].y += (dv * dH1X.y);
    water_a[wi].z -= (dv * dH1X.z);
    water_b[wj].z += (dv * dH1X.z);

    // O - H2
    r2 = pow(rH2X, 2);
    Vel = Coul * crg_ow * crg_hw * rH2X;
    *Ecoul += Vel;
    dv = r2 * (-Vel);
    j += 1; //point to H2 in j-molecule
    wj += 1;

    water_a[wi].x -= (dv * dH2X.x);
    water_b[wj].x += (dv * dH2X.x);
    water_a[wi].y -= (dv * dH2X.y);
    water_b[wj].y += (dv * dH2X.y);
    water_a[wi].z -= (dv * dH2X.z);
    water_b[wj].z += (dv * dH2X.z);

    // --- H1 - (O,H1,H2) ---
    i += 1; //Point to H1 in i-molecule
    wi += 1;
    j -= 2; //Point to O in j-molecule
    wj -= 2;

    // H1 - O (X=O)
    dOX.x = Ys[j].x - Xs[i].x;
    dOX.y = Ys[j].y - Xs[i].y;
    dOX.z = Ys[j].z - Xs[i].z;
    rOX = pow(dOX.x, 2) + pow(dOX.y, 2) + pow(dOX.z, 2);

    // H1 - H1 (X=H1)
    j += 1;
    wj += 1;

    dH1X.x = Ys[j].x - Xs[i].x;
    dH1X.y = Ys[j].y - Xs[i].y;
    dH1X.z = Ys[j].z - Xs[i].z;
    rH1X = pow(dH1X.x, 2) + pow(dH1X.y, 2) + pow(dH1X.z, 2);

    // H1 - H2 (X=H2)
    j += 1;
    wj += 1;

    dH2X.x = Ys[j].x - Xs[i].x;
    dH2X.y = Ys[j].y - Xs[i].y;
    dH2X.z = Ys[j].z - Xs[i].z;
    rH2X = pow(dH2X.x, 2) + pow(dH2X.y, 2) + pow(dH2X.z, 2);
    rOX = sqrt(1 / rOX);
    rH1X = sqrt(1 / rH1X);
    rH2X = sqrt(1 / rH2X);

    // H1 - O
    r2 = rOX * rOX;
    Vel = Coul * crg_hw * crg_ow * rOX;
    *Ecoul += Vel;
    dv = r2 * (-Vel);
    j -= 2; //move pointer back to O in interacting molecule
    wj -= 2;
    water_a[wi].x -= (dv * dOX.x);
    water_b[wj].x += (dv * dOX.x);
    water_a[wi].y -= (dv * dOX.y);
    water_b[wj].y += (dv * dOX.y);
    water_a[wi].z -= (dv * dOX.z);
    water_b[wj].z += (dv * dOX.z);

    // H1 - H1
    r2 = pow(rH1X, 2);
    Vel = Coul * crg_hw * crg_hw * rH1X;
    *Ecoul += Vel;
    dv = r2 * (-Vel);
    j += 1; //point to H1 in j-molecule
    wj += 1;
    water_a[wi].x -= (dv * dH1X.x);
    water_b[wj].x += (dv * dH1X.x);
    water_a[wi].y -= (dv * dH1X.y);
    water_b[wj].y += (dv * dH1X.y);
    water_a[wi].z -= (dv * dH1X.z);
    water_b[wj].z += (dv * dH1X.z);

    // H1 - H2
    r2 = pow(rH2X, 2);
    Vel = Coul * crg_hw * crg_hw * rH2X;
    *Ecoul += Vel;
    dv = r2 * (-Vel);
    j += 1; //point to H2 in j-molecule
    wj += 1;
    water_a[wi].x -= (dv * dH2X.x);
    water_b[wj].x += (dv * dH2X.x);
    water_a[wi].y -= (dv * dH2X.y);
    water_b[wj].y += (dv * dH2X.y);
    water_a[wi].z -= (dv * dH2X.z);
    water_b[wj].z += (dv * dH2X.z);

    // --- H2 - (O,H1,H2) ---
    i += 1; //Point to H2 in i-molecule
    wi += 1;
    j -= 2; //Point to O in j-molecule
    wj -= 2;

    // H2 - O (X=O)
    dOX.x = Ys[j].x - Xs[i].x;
    dOX.y = Ys[j].y - Xs[i].y;
    dOX.z = Ys[j].z - Xs[i].z;
    rOX = pow(dOX.x, 2) + pow(dOX.y, 2) + pow(dOX.z, 2);

    // H2 - H1 (X=H1)
    j += 1;
    wj += 1;

    dH1X.x = Ys[j].x - Xs[i].x;
    dH1X.y = Ys[j].y - Xs[i].y;
    dH1X.z = Ys[j].z - Xs[i].z;
    rH1X = pow(dH1X.x, 2) + pow(dH1X.y, 2) + pow(dH1X.z, 2);

    // H2 - H2 (X=H2)
    j += 1;
    wj += 1;

    dH2X.x = Ys[j].x - Xs[i].x;
    dH2X.y = Ys[j].y - Xs[i].y;
    dH2X.z = Ys[j].z - Xs[i].z;
    rH2X = pow(dH2X.x, 2) + pow(dH2X.y, 2) + pow(dH2X.z, 2);
    rOX = sqrt(1 / rOX);
    rH1X = sqrt(1 / rH1X);
    rH2X = sqrt(1 / rH2X);

    // H2 - O
    r2 = rOX * rOX;
    Vel = Coul * crg_hw * crg_ow * rOX;
    *Ecoul += Vel;
    dv = r2 * (-Vel);
    j -= 2; //move pointer back to O in interacting molecule
    wj -= 2;
    water_a[wi].x -= (dv * dOX.x);
    water_b[wj].x += (dv * dOX.x);
    water_a[wi].y -= (dv * dOX.y);
    water_b[wj].y += (dv * dOX.y);
    water_a[wi].z -= (dv * dOX.z);
    water_b[wj].z += (dv * dOX.z);

    // H2 - H1
    r2 = pow(rH1X, 2);
    Vel = Coul * crg_hw * crg_hw * rH1X;
    *Ecoul += Vel;
    dv = r2 * (-Vel);
    j += 1; //point to H1 in j-molecule
    wj += 1;
    water_a[wi].x -= (dv * dH1X.x);
    water_b[wj].x += (dv * dH1X.x);
    water_a[wi].y -= (dv * dH1X.y);
    water_b[wj].y += (dv * dH1X.y);
    water_a[wi].z -= (dv * dH1X.z);
    water_b[wj].z += (dv * dH1X.z);

    // H1 - H2
    r2 = pow(rH2X, 2);
    Vel = Coul * crg_hw * crg_hw * rH2X;
    *Ecoul += Vel;
    dv = r2 * (-Vel);
    j += 1; //point to H2 in j-molecule
    wj += 1;
    water_a[wi].x -= (dv * dH2X.x);
    water_b[wj].x += (dv * dH2X.x);
    water_a[wi].y -= (dv * dH2X.y);
    water_b[wj].y += (dv * dH2X.y);
    water_a[wi].z -= (dv * dH2X.z);
    water_b[wj].z += (dv * dH2X.z);
}

__global__ void calc_ww_dvel_matrix(int n_waters, double crg_ow, double crg_hw, double A_OO, double B_OO,
    coord_t *X, double *Evdw, double *Ecoul, calc_ww_t *MAT) {
    // Block index
    int bx = blockIdx.x;
    int by = blockIdx.y;

    // Thread index
    int tx = threadIdx.x;
    int ty = threadIdx.y;

    // if (bx == 9 && by == 0) printf("bx = %d by = %d tx = %d ty = %d\n", bx, by, tx, ty);

    int aStart = 3 * BLOCK_SIZE * by;
    int bStart = 3 * BLOCK_SIZE * bx;

    if (aStart + 3 * ty >= 3 * n_waters) return;
    if (bStart + 3 * tx >= 3 * n_waters) return;

    __shared__ coord_t Xs[3 * BLOCK_SIZE];
    __shared__ coord_t Ys[3 * BLOCK_SIZE];
    
    if (tx == 0) {
        Xs[3 * ty    ] = X[aStart + 3 * ty    ];
        Xs[3 * ty + 1] = X[aStart + 3 * ty + 1];
        Xs[3 * ty + 2] = X[aStart + 3 * ty + 2];
    }

    if (ty == 0) {
        Ys[3 * tx    ] = X[bStart + 3 * tx    ];
        Ys[3 * tx + 1] = X[bStart + 3 * tx + 1];
        Ys[3 * tx + 2] = X[bStart + 3 * tx + 2];
    }

    __syncthreads();

    if (bx < by || (bx == by && tx < ty)) return;

    dvel_t water_a[3], water_b[3];
    memset(&water_a, 0, 3 * sizeof(dvel_t));
    memset(&water_b, 0, 3 * sizeof(dvel_t));

    int row = by * BLOCK_SIZE + ty;
    int column = bx * BLOCK_SIZE + tx;
    
    // if (row == 0 && column == 1) {
    //     printf("Xs[0] = %f\n", Xs[0]);
    //     printf("Ys[0] = %f\n", Ys[0]);
    //     printf("Xs[1] = %f\n", Xs[1]);
    //     printf("Ys[1] = %f\n", Ys[1]);
    //     printf("Xs[2] = %f\n", Xs[2]);
    //     printf("Ys[2] = %f\n", Ys[2]);
    //     printf("Xs[3] = %f\n", Xs[3]);
    //     printf("Ys[3] = %f\n", Ys[3]);
    //     printf("Xs[4] = %f\n", Xs[4]);
    //     printf("Ys[4] = %f\n", Ys[4]);

    //     printf("Ys[%d] = %f Xs[%d] = %f\n", 3 * ty, Ys[3 * ty], 3 * tx, Xs[3 * tx]);
    // }

    if (bx != by || tx != ty) {
        double evdw, ecoul;
        calc_ww_dvel_matrix_incr(ty, tx, crg_ow, crg_hw, A_OO, B_OO, Xs, Ys, &evdw, &ecoul, water_a, water_b);
    }

    // if (row == 0 && column == 1) {
    //     printf("water_a = %f %f %f water_b = %f %f %f\n", water_a[0].x, water_a[0].y, water_a[0].z, water_b[0].x, water_b[0].y, water_b[0].z);
    // }

    set_water(n_waters, row, column, water_a, MAT);
    set_water(n_waters, column, row, water_b, MAT);

    __syncthreads();
}

__global__ void calc_ww_dvel_vector(int n_waters, dvel_t *DV, calc_ww_t *MAT) {
    int row = blockIdx.x*blockDim.x + threadIdx.x;
    if (row >= n_waters) return;

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

    for (int i = 0; i < n_waters; i++) {
        if (i != row) {
            dO.x += MAT[i + n_waters * row].O.x;
            dO.y += MAT[i + n_waters * row].O.y;
            dO.z += MAT[i + n_waters * row].O.z;
            dH1.x += MAT[i + n_waters * row].H1.x;
            dH1.y += MAT[i + n_waters * row].H1.y;
            dH1.z += MAT[i + n_waters * row].H1.z;
            dH2.x += MAT[i + n_waters * row].H2.x;
            dH2.y += MAT[i + n_waters * row].H2.y;
            dH2.z += MAT[i + n_waters * row].H2.z;
        }
    }

    DV[3*row].x += dO.x;
    DV[3*row].y += dO.y;
    DV[3*row].z += dO.z;
    DV[3*row+1].x += dH1.x;
    DV[3*row+1].y += dH1.y;
    DV[3*row+1].z += dH1.z;
    DV[3*row+2].x += dH2.x;
    DV[3*row+2].y += dH2.y;
    DV[3*row+2].z += dH2.z;

    __syncthreads();
}

__device__ void calc_pw_dvel_matrix_incr(int row, int pi, int column, int j, int n_patoms,
    coord_t *Ps, coord_t *Xs, double *Evdw, double *Ecoul, calc_pw_t *pw,
    ccharge_t *D_ccharges, charge_t *D_charges, catype_t *D_catypes, atype_t *D_atypes, p_atom_t *D_patoms) {

    coord_t da;
    double r2a, ra, r6a;
    double Vela, V_a, V_b;
    double dva;
    double qi, qj;
    double ai_aii, aj_aii, ai_bii, aj_bii;
    catype_t ai_type, aj_type;
    int i;
    atype_t i_type, j_type;

    i = D_patoms[pi].a-1;
    // if (excluded[i] || excluded[j]) continue;
    qi = D_ccharges[D_charges[i].code - 1].charge;
    qj = D_ccharges[D_charges[n_patoms + j].code - 1].charge; //TODO: FIX THIS!!! WILL NOT WORK WITH QATOMS!!!!!

    // if (pi < 100 && j < 100){
    //     printf("qi = %f qj = %f\n", qi, qj);
    // }

    ai_type = D_catypes[D_atypes[i].code - 1];
    aj_type = D_catypes[D_atypes[n_patoms + j].code - 1];

    da.x = Xs[column].x - Ps[row].x;
    da.y = Xs[column].y - Ps[row].y;
    da.z = Xs[column].z - Ps[row].z;
    r2a = 1 / (pow(da.x, 2) + pow(da.y, 2) + pow(da.z, 2));
    ra = sqrt(r2a);
    r6a = r2a * r2a * r2a;

    Vela = Coul * qi * qj * ra;

    ai_aii = ai_type.aii_normal;
    aj_aii = aj_type.aii_normal;
    ai_bii = ai_type.bii_normal;
    aj_bii = aj_type.bii_normal;

    V_a = r6a * r6a * ai_aii * aj_aii;
    V_b = r6a * ai_bii * aj_bii;
    dva = r2a * ( -Vela -12 * V_a + 6 * V_b);


    pw->P.x -= dva * da.x;
    pw->P.y -= dva * da.y;
    pw->P.z -= dva * da.z;

    pw->W.x += dva * da.x;
    pw->W.y += dva * da.y;
    pw->W.z += dva * da.z;

    *Evdw += Vela;
    *Ecoul += (V_a - V_b);

    if (pi == 522 && j == 175) {
        printf("Vela = %f V_a = %f V_b = %f P = %f %f %f ai_aii = %f aj_aii = %f\n", Vela, V_a, V_b, pw->P.x, pw->P.y, pw->P.z, ai_aii, aj_aii);
    }

    // if (pi < 100 && j < 100) printf("Evdw = %f Ecoul = %f\n", *Evdw, *Ecoul);
}

__global__ void calc_pw_dvel_matrix(int n_patoms, int n_waters,
    coord_t *P, coord_t *X, double *Evdw, double *Ecoul, calc_pw_t *PW_MAT,
    ccharge_t *D_ccharges, charge_t *D_charges, catype_t *D_catypes, atype_t *D_atypes, p_atom_t *D_patoms) {
    // Block index
    int bx = blockIdx.x;
    int by = blockIdx.y;

    // Thread index
    int tx = threadIdx.x;
    int ty = threadIdx.y;

    // if (bx == 0 && by == 0) printf("bx = %d by = %d tx = %d ty = %d\n", bx, by, tx, ty);

    int aStart = BLOCK_SIZE * by;
    int bStart = BLOCK_SIZE * bx;

    if (aStart + ty >= n_patoms) return;
    if (bStart + tx >= 3 * n_waters) return;

    // if (bx == 8 && by == 1) printf("bx = %d by = %d tx = %d ty = %d\n", bx, by, tx, ty);

    __shared__ coord_t Ps[BLOCK_SIZE];
    __shared__ coord_t Xs[BLOCK_SIZE];

    Ps[ty] = P[D_patoms[aStart + ty].a-1];
    Xs[tx] = X[bStart + tx];

    __syncthreads();

    calc_pw_t pw;
    memset(&pw, 0, sizeof(calc_pw_t));

    int row = by * BLOCK_SIZE + ty;
    int column = bx * BLOCK_SIZE + tx;

    // if (row == 0 && column == 1) {
    //     printf("Xs[0] = %f\n", Xs[0]);
    //     printf("Ys[0] = %f\n", Ys[0]);
    //     printf("Xs[1] = %f\n", Xs[1]);
    //     printf("Ys[1] = %f\n", Ys[1]);
    //     printf("Xs[2] = %f\n", Xs[2]);
    //     printf("Ys[2] = %f\n", Ys[2]);
    //     printf("Xs[3] = %f\n", Xs[3]);
    //     printf("Ys[3] = %f\n", Ys[3]);
    //     printf("Xs[4] = %f\n", Xs[4]);
    //     printf("Ys[4] = %f\n", Ys[4]);

    //     printf("Ys[%d] = %f Xs[%d] = %f\n", 3 * ty, Ys[3 * ty], 3 * tx, Xs[3 * tx]);
    // }

    // if (bx == 8 && by == 1) printf("bx = %d by = %d tx = %d ty = %d\n", bx, by, tx, ty);
    // __device__ void calc_pw_dvel_matrix_incr(int row, int pi, int column, int j, int n_patoms,
    // coord_t *Ps, coord_t *Xs, double *Evdw, double *Ecoul, calc_pw_t *pw,
    // ccharge_t *D_ccharges, charge_t *D_charges, catype_t *D_catypes, atype_t *D_atypes, p_atom_t *D_patoms)
    double evdw, ecoul;
    calc_pw_dvel_matrix_incr(ty, aStart + ty, tx, bStart + tx, n_patoms, Ps, Xs, &evdw, &ecoul, &pw, D_ccharges, D_charges, D_catypes, D_atypes, D_patoms);

    // if (row == 0 && column == 1) {
    //     printf("water_a = %f %f %f water_b = %f %f %f\n", water_a[0].x, water_a[0].y, water_a[0].z, water_b[0].x, water_b[0].y, water_b[0].z);
    // }

    // if (bx == 8 && by == 1) printf("n_qatoms = %d\n", n_qatoms);
    // if (bx == 8 && by == 1) printf("qi = %d j = %d charge[%d] = %f\n", row, column, row + n_qatoms, D_qcharges[row + n_qatoms * 1].q);

    PW_MAT[column + 3 * n_waters * row] = pw;

    __syncthreads();
}

__global__ void calc_pw_dvel_vector_row(int n_patoms, int n_waters, dvel_t *DV_P, dvel_t *DV, calc_pw_t *PW_MAT) {
    int row = blockIdx.x*blockDim.x + threadIdx.x;
    if (row >= n_patoms) return;

    dvel_t dP;

    dP.x = 0;
    dP.y = 0;
    dP.z = 0;

    for (int i = 0; i < 3*n_waters; i++) {
        dP.x += PW_MAT[i + 3*n_waters * row].P.x;
        dP.y += PW_MAT[i + 3*n_waters * row].P.y;
        dP.z += PW_MAT[i + 3*n_waters * row].P.z;
    }

    DV_P[row].x += dP.x;
    DV_P[row].y += dP.y;
    DV_P[row].z += dP.z;

    __syncthreads();
}

__global__ void calc_pw_dvel_vector_column(int n_patoms, int n_waters, dvel_t *DV_P, dvel_t *DV, calc_pw_t *PW_MAT) {
    int column = blockIdx.x*blockDim.x + threadIdx.x;
    if (column >= 3*n_waters) return;

    dvel_t dW;

    dW.x = 0;
    dW.y = 0;
    dW.z = 0;

    for (int i = 0; i < n_patoms; i++) {
        dW.x += PW_MAT[column + 3*n_waters * i].W.x;
        dW.y += PW_MAT[column + 3*n_waters * i].W.y;
        dW.z += PW_MAT[column + 3*n_waters * i].W.z;
        // if (column == 1754) printf("dvelocities[%d].x += %f (MAX %d)\n", column + n_patoms, dW.x, n_patoms + n_waters*3 - 1);
    }


    DV[column].x += dW.x;
    DV[column].y += dW.y;
    DV[column].z += dW.z;

    __syncthreads();
}

void clean_solvent() {
    cudaFree(P);
    cudaFree(X);
    cudaFree(WW_MAT);
    cudaFree(PW_MAT);
    cudaFree(DV);
    cudaFree(DV_P);

    cudaFree(D_ccharges);
    cudaFree(D_charges);
    cudaFree(D_catypes);
    cudaFree(D_atypes);
    cudaFree(D_patoms);

    free(h_WW_MAT);
    free(h_PW_MAT);

    printf("All memory freed\n");
}