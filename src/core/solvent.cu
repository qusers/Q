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
dvel_t *MAT, *DV;

//ONLY call if there are actually solvent atoms, or get segfaulted
void calc_nonbonded_ww_forces() {
    double rOX, rH1X, rH2X, r2;
    coord_t dOX, dH1X, dH2X;
    double Vel, V_a, V_b, dv;

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
}

void calc_nonbonded_ww_forces_host() {
    int mem_size_X = 3 * n_waters * sizeof(coord_t);
    int mem_size_MAT = 3 * 3 * n_waters * n_waters * sizeof(dvel_t);
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

        cudaMalloc((void**) &X, mem_size_X);
        error = cudaMalloc((void**) &MAT, mem_size_MAT);
        if (error != cudaSuccess) {
            printf(">>> FATAL: memory for matrix could not be allocated. Exiting...\n");
            exit(EXIT_FAILURE);
        }
        cudaMalloc((void**) &DV, mem_size_DV);    
    }

    clock_t start = clock();

    cudaMemcpy(X, coords, mem_size_X, cudaMemcpyHostToDevice);
    cudaMemset(MAT, 0, mem_size_MAT);
    cudaMemcpy(DV, dvelocities, mem_size_DV, cudaMemcpyHostToDevice);

    clock_t end = clock();

    printf("Time to copy data over: %f\n", (end-start) / (double)CLOCKS_PER_SEC);

    dim3 threads,grid;

    threads = dim3(BLOCK_SIZE, BLOCK_SIZE);
    grid = dim3(3 * n_waters / threads.x, 3 * n_waters / threads.y);

    double evdw, ecoul;

    calc_ww_dvel_matrix<<<grid, threads>>>(n_waters, crg_ow, crg_hw, A_OO, B_OO, X, &evdw, &ecoul, MAT);
    calc_ww_dvel_vector_rows<<<(3 * 3 * n_waters * n_waters / BLOCK_SIZE), BLOCK_SIZE>>>(n_waters, DV, MAT);
    calc_ww_dvel_vector_columns<<<(3 * 3 * n_waters * n_waters / BLOCK_SIZE), BLOCK_SIZE>>>(n_waters, DV, MAT);

    cudaMemcpy(dvelocities, DV, mem_size_DV, cudaMemcpyDeviceToHost);
}

__device__ void calc_ww_dvel_matrix_incr(int row, int column, double crg_ow, double crg_hw, double A_OO, double B_OO,
    coord_t *Xs, coord_t *Ys, double *Evdw, double *Ecoul, dvel_t *water_a, dvel_t *water_b) {
    
    // Do stuff
    double rOX, rH1X, rH2X, r2;
    coord_t dOX, dH1X, dH2X;
    double Vel, V_a, V_b, dv;

    int w_i = 0;
    int w_j = 0;

    int i = 3 * row;
    int j = 3 * column;

    dOX.x = Ys[j].x - Xs[i].x;
    dOX.y = Ys[j].y - Xs[i].y;
    dOX.z = Ys[j].z - Xs[i].z;
    rOX = pow(dOX.x, 2) + pow(dOX.y, 2) + pow(dOX.z, 2);

    // O - H1
    j += 1;
    w_j++;

    dH1X.x = Ys[j].x - Xs[i].x;
    dH1X.y = Ys[j].y - Xs[i].y;
    dH1X.z = Ys[j].z - Xs[i].z;
    rH1X = pow(dH1X.x, 2) + pow(dH1X.y, 2) + pow(dH1X.z, 2);

    // O-H2 (X=H2)
    j += 1;
    w_j++;

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
    w_j -= 2;
    water_a[w_i].x -= (dv * dOX.x);
    water_b[w_j].x += (dv * dOX.x);
    water_a[w_i].y -= (dv * dOX.y);
    water_b[w_j].y += (dv * dOX.y);
    water_a[w_i].z -= (dv * dOX.z);
    water_b[w_j].z += (dv * dOX.z);

    // O - H1
    r2 = pow(rH1X, 2);
    Vel = Coul * crg_ow * crg_hw * rH1X;
    *Ecoul += Vel;
    dv = r2 * (-Vel);
    j += 1; //point to H1 in j-molecule
    w_j++;
    water_a[w_i].x -= (dv * dH1X.x);
    water_b[w_j].x += (dv * dH1X.x);
    water_a[w_i].y -= (dv * dH1X.y);
    water_b[w_j].y += (dv * dH1X.y);
    water_a[w_i].z -= (dv * dH1X.z);
    water_b[w_j].z += (dv * dH1X.z);

    // O - H2
    r2 = pow(rH2X, 2);
    Vel = Coul * crg_ow * crg_hw * rH2X;
    *Ecoul += Vel;
    dv = r2 * (-Vel);
    j += 1; //point to H2 in j-molecule
    w_j++;
    water_a[w_i].x -= (dv * dH2X.x);
    water_b[w_j].x += (dv * dH2X.x);
    water_a[w_i].y -= (dv * dH2X.y);
    water_b[w_j].y += (dv * dH2X.y);
    water_a[w_i].z -= (dv * dH2X.z);
    water_b[w_j].z += (dv * dH2X.z);

    // --- H1 - (O,H1,H2) ---
    i += 1; //Point to H1 in i-molecule
    w_i++;
    j -= 2; //Point to O in j-molecule
    w_j -= 2;

    // H1 - O (X=O)
    dOX.x = Ys[j].x - Xs[i].x;
    dOX.y = Ys[j].y - Xs[i].y;
    dOX.z = Ys[j].z - Xs[i].z;
    rOX = pow(dOX.x, 2) + pow(dOX.y, 2) + pow(dOX.z, 2);

    // H1 - H1 (X=H1)
    j += 1;
    w_j++;

    dH1X.x = Ys[j].x - Xs[i].x;
    dH1X.y = Ys[j].y - Xs[i].y;
    dH1X.z = Ys[j].z - Xs[i].z;
    rH1X = pow(dH1X.x, 2) + pow(dH1X.y, 2) + pow(dH1X.z, 2);

    // H1 - H2 (X=H2)
    j += 1;
    w_j++;

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
    w_j -= 2;
    water_a[w_i].x -= (dv * dOX.x);
    water_b[w_j].x += (dv * dOX.x);
    water_a[w_i].y -= (dv * dOX.y);
    water_b[w_j].y += (dv * dOX.y);
    water_a[w_i].z -= (dv * dOX.z);
    water_b[w_j].z += (dv * dOX.z);

    // H1 - H1
    r2 = pow(rH1X, 2);
    Vel = Coul * crg_hw * crg_hw * rH1X;
    *Ecoul += Vel;
    dv = r2 * (-Vel);
    j += 1; //point to H1 in j-molecule
    w_j++;
    water_a[w_i].x -= (dv * dH1X.x);
    water_b[w_j].x += (dv * dH1X.x);
    water_a[w_i].y -= (dv * dH1X.y);
    water_b[w_j].y += (dv * dH1X.y);
    water_a[w_i].z -= (dv * dH1X.z);
    water_b[w_j].z += (dv * dH1X.z);

    // H1 - H2
    r2 = pow(rH2X, 2);
    Vel = Coul * crg_hw * crg_hw * rH2X;
    *Ecoul += Vel;
    dv = r2 * (-Vel);
    j += 1; //point to H2 in j-molecule
    w_j++;
    water_a[w_i].x -= (dv * dH2X.x);
    water_b[w_j].x += (dv * dH2X.x);
    water_a[w_i].y -= (dv * dH2X.y);
    water_b[w_j].y += (dv * dH2X.y);
    water_a[w_i].z -= (dv * dH2X.z);
    water_b[w_j].z += (dv * dH2X.z);

    // --- H2 - (O,H1,H2) ---
    i += 1; //Point to H2 in i-molecule
    w_i++;
    j -= 2; //Point to O in j-molecule
    w_j -= 2;

    // H2 - O (X=O)
    dOX.x = Ys[j].x - Xs[i].x;
    dOX.y = Ys[j].y - Xs[i].y;
    dOX.z = Ys[j].z - Xs[i].z;
    rOX = pow(dOX.x, 2) + pow(dOX.y, 2) + pow(dOX.z, 2);

    // H2 - H1 (X=H1)
    j += 1;
    w_j++;

    dH1X.x = Ys[j].x - Xs[i].x;
    dH1X.y = Ys[j].y - Xs[i].y;
    dH1X.z = Ys[j].z - Xs[i].z;
    rH1X = pow(dH1X.x, 2) + pow(dH1X.y, 2) + pow(dH1X.z, 2);

    // H2 - H2 (X=H2)
    j += 1;
    w_j++;

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
    w_j -= 2;
    water_a[w_i].x -= (dv * dOX.x);
    water_b[w_j].x += (dv * dOX.x);
    water_a[w_i].y -= (dv * dOX.y);
    water_b[w_j].y += (dv * dOX.y);
    water_a[w_i].z -= (dv * dOX.z);
    water_b[w_j].z += (dv * dOX.z);

    // H2 - H1
    r2 = pow(rH1X, 2);
    Vel = Coul * crg_hw * crg_hw * rH1X;
    *Ecoul += Vel;
    dv = r2 * (-Vel);
    j += 1; //point to H1 in j-molecule
    w_j++;
    water_a[w_i].x -= (dv * dH1X.x);
    water_b[w_j].x += (dv * dH1X.x);
    water_a[w_i].y -= (dv * dH1X.y);
    water_b[w_j].y += (dv * dH1X.y);
    water_a[w_i].z -= (dv * dH1X.z);
    water_b[w_j].z += (dv * dH1X.z);

    // H1 - H2
    r2 = pow(rH2X, 2);
    Vel = Coul * crg_hw * crg_hw * rH2X;
    *Ecoul += Vel;
    dv = r2 * (-Vel);
    j += 1; //point to H2 in j-molecule
    w_j++;
    water_a[w_i].x -= (dv * dH2X.x);
    water_b[w_j].x += (dv * dH2X.x);
    water_a[w_i].y -= (dv * dH2X.y);
    water_b[w_j].y += (dv * dH2X.y);
    water_a[w_i].z -= (dv * dH2X.z);
    water_b[w_j].z += (dv * dH2X.z);
}

__global__ void calc_ww_dvel_matrix(int n_waters, double crg_ow, double crg_hw, double A_OO, double B_OO,
    coord_t *X, double *Evdw, double *Ecoul, dvel_t *MAT) {
    // Block index
    int bx = blockIdx.x;
    int by = blockIdx.y;

    // Thread index
    int tx = threadIdx.x;
    int ty = threadIdx.y;
    
    __shared__ coord_t Xs[3 * BLOCK_SIZE];
    __shared__ coord_t Ys[3 * BLOCK_SIZE];

    int aStart = 3 * BLOCK_SIZE * by;
    int bStart = 3 * BLOCK_SIZE * bx;

    if (aStart + 3 * tx >= 3 * n_waters) return;
    if (bStart + 3 * ty >= 3 * n_waters) return;

    if (bx == by) {
        Xs[3 * tx    ] = X[aStart + 3 * tx    ];
        Xs[3 * tx + 1] = X[aStart + 3 * tx + 1];
        Xs[3 * tx + 2] = X[aStart + 3 * tx + 2];
    }
    else {
        Xs[3 * tx    ] = X[aStart + 3 * tx    ];
        Xs[3 * tx + 1] = X[aStart + 3 * tx + 1];
        Xs[3 * tx + 2] = X[aStart + 3 * tx + 2];
    
        Ys[3 * ty    ] = X[bStart + 3 * ty    ];
        Ys[3 * ty + 1] = X[bStart + 3 * ty + 1];
        Ys[3 * ty + 2] = X[bStart + 3 * ty + 2];
    }

    __syncthreads();

    if (bx == by && tx == ty) {
        Ys[3 * ty    ] = Xs[3 * ty    ];
        Ys[3 * ty + 1] = Xs[3 * ty + 1];
        Ys[3 * ty + 2] = Xs[3 * ty + 2];

    }

    __syncthreads();

    dvel_t water_a[3], water_b[3];

    if (bx != by || tx != ty) {
        double evdw, ecoul;
        calc_ww_dvel_matrix_incr(tx, ty, crg_ow, crg_hw, A_OO, B_OO, Xs, Ys, &evdw, &ecoul, water_a, water_b);
    }

    __syncthreads();

    int ix = n_waters * 3 * BLOCK_SIZE * by + BLOCK_SIZE * bx;

    MAT[ix + n_waters * 3 * ty + tx    ] = water_a[0];
    MAT[ix + n_waters * 3 * ty + tx + 1] = water_a[1];
    MAT[ix + n_waters * 3 * ty + tx + 2] = water_a[2];

    __syncthreads();
}

__global__ void calc_ww_dvel_vector_rows(int n_waters, dvel_t *DV, dvel_t *MAT) {
    int row = blockIdx.x*blockDim.x + threadIdx.x;
    if (row >= 3 * n_waters) return;

    dvel_t dv;
    dv.x = 0;
    dv.y = 0;
    dv.z = 0;

    for (int i = 0; i < 3 * n_waters; i++) {
        if (i > row) {
            dv.x += MAT[row + 3 * n_waters * i].x;
            dv.y += MAT[row + 3 * n_waters * i].y;
            dv.z += MAT[row + 3 * n_waters * i].z;
        }
    }

    DV[row].x += dv.x;
    DV[row].y += dv.y;
    DV[row].z += dv.z;

    __syncthreads();
}

__global__ void calc_ww_dvel_vector_columns(int n_waters, dvel_t *DV, dvel_t *MAT) {
    int column = blockIdx.x*blockDim.x + threadIdx.x;
    if (column >= 3 * n_waters) return;

    dvel_t dv;
    dv.x = 0;
    dv.y = 0;
    dv.z = 0;

    for (int i = 0; i < 3 * n_waters; i++) {
        if (i < column) {
            dv.x += MAT[i + 3 * n_waters * column].x;
            dv.y += MAT[i + 3 * n_waters * column].y;
            dv.z += MAT[i + 3 * n_waters * column].z;
        }
    }

    DV[column].x -= dv.x;
    DV[column].y -= dv.y;
    DV[column].z -= dv.z;

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