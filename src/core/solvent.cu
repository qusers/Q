#include "system.h"
#include "solvent.h"

#include <stdio.h>
#include <time.h>
#include <unistd.h>

/* =============================================
 * == SOLVENT INTERACTIONS
 * =============================================
 */

coord_t *X;
dvel_t *DV, *h_DV;
calc_water_t *MAT, *h_MAT;

//ONLY call if there are actually solvent atoms, or get segfaulted
void calc_nonbonded_ww_forces() {
    energies.Ucoul = 0;
    energies.Uvdw = 0;

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
            dvel_t water_a[3], water_b[3];
            memset(&water_a, 0, 3 * sizeof(dvel_t));
            memset(&water_b, 0, 3 * sizeof(dvel_t));
                
            calc_nonbonded_ww_forces_incr(i/3, j/3, crg_ow, crg_hw, A_OO, B_OO, coords, coords, &energies.Uvdw, &energies.Ucoul, water_a, water_b);

            // printf("MAT[%d][%d].O = %f %f %f\n", i/3, j/3, water_a[0].x, water_a[0].y, water_a[0].z);
            // printf("MAT[%d][%d].O = %f %f %f\n", j/3, i/3, water_b[0].x, water_b[0].y, water_b[0].z);

            dvelocities[i].x += water_a[0].x;
            dvelocities[i].y += water_a[0].y;
            dvelocities[i].z += water_a[0].z;
            dvelocities[i+1].x += water_a[1].x;
            dvelocities[i+1].y += water_a[1].y;
            dvelocities[i+1].z += water_a[1].z;
            dvelocities[i+2].x += water_a[2].x;
            dvelocities[i+2].y += water_a[2].y;
            dvelocities[i+2].z += water_a[2].z;

            dvelocities[j].x += water_b[0].x;
            dvelocities[j].y += water_b[0].y;
            dvelocities[j].z += water_b[0].z;
            dvelocities[j+1].x += water_b[1].x;
            dvelocities[j+1].y += water_b[1].y;
            dvelocities[j+1].z += water_b[1].z;
            dvelocities[j+2].x += water_b[2].x;
            dvelocities[j+2].y += water_b[2].y;
            dvelocities[j+2].z += water_b[2].z;
        }
    }
}

void calc_nonbonded_ww_forces_incr(int row, int column, double crg_ow, double crg_hw, double A_OO, double B_OO,
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

void calc_nonbonded_ww_forces_host() {
    int mem_size_X = 3 * n_waters * sizeof(coord_t);
    int mem_size_MAT = n_waters * n_waters * sizeof(calc_water_t);
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

        cudaSetDevice(1);

        cudaMalloc((void**) &X, mem_size_X);
        error = cudaMalloc((void**) &MAT, mem_size_MAT);
        if (error != cudaSuccess) {
            printf(">>> FATAL: memory for matrix could not be allocated. Exiting...\n");
            exit(EXIT_FAILURE);
        }
        cudaMalloc((void**) &DV, mem_size_DV);    
        h_DV = (dvel_t*) malloc(3 * n_waters * sizeof(dvel_t));
        h_MAT = (calc_water_t*) malloc(mem_size_MAT);
    }

    cudaMemcpy(X, coords, mem_size_X, cudaMemcpyHostToDevice);
    cudaMemcpy(DV, dvelocities, mem_size_DV, cudaMemcpyHostToDevice);

    dim3 threads,grid;

    threads = dim3(BLOCK_SIZE, BLOCK_SIZE);
    grid = dim3((n_waters + BLOCK_SIZE - 1) / threads.x, (n_waters + BLOCK_SIZE - 1) / threads.y);

    double evdw, ecoul;

    calc_ww_dvel_matrix<<<grid, threads>>>(n_waters, crg_ow, crg_hw, A_OO, B_OO, X, &evdw, &ecoul, MAT);
    calc_ww_dvel_vector_rows<<<((n_waters+BLOCK_SIZE - 1) / BLOCK_SIZE), BLOCK_SIZE>>>(n_waters, DV, MAT);

    // cudaMemcpy(h_MAT, MAT, mem_size_MAT, cudaMemcpyDeviceToHost);

    // for (int i = 0; i < n_waters; i++) {
    //     for (int j = 0; j < n_waters; j++) {
    //         printf("MAT[%d][%d].O = %f %f %f\n", i, j, h_MAT[i * n_waters + j].O.x, h_MAT[i * n_waters + j].O.y, h_MAT[i * n_waters + j].O.z);
    //     }
    // }

    cudaMemcpy(dvelocities, DV, mem_size_DV, cudaMemcpyDeviceToHost);
}

__device__ void set_water(int n_waters, int row, int column, dvel_t *val, calc_water_t *MAT) {
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
    coord_t *X, double *Evdw, double *Ecoul, calc_water_t *MAT) {
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

    if (bx != by || tx != ty) {
        double evdw, ecoul;
        calc_ww_dvel_matrix_incr(ty, tx, crg_ow, crg_hw, A_OO, B_OO, Xs, Ys, &evdw, &ecoul, water_a, water_b);
    }

    int row = by * BLOCK_SIZE + ty;
    int column = bx * BLOCK_SIZE + tx;

    set_water(n_waters, row, column, water_a, MAT);
    set_water(n_waters, column, row, water_b, MAT);

    __syncthreads();
}

__global__ void calc_ww_dvel_vector_rows(int n_waters, dvel_t *DV, calc_water_t *MAT) {
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
        }
    }
}

void clean_solvent() {
    cudaFree(X);
    cudaFree(MAT);
    cudaFree(DV);
}
