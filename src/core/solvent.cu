#include "system.h"
#include "solvent.h"
#include <stdio.h>

/* =============================================
 * == SOLVENT INTERACTIONS
 * =============================================
 */

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
            E_nonbond_ww.Uvdw += (V_a - V_b);
            E_nonbond_ww.Ucoul += Vel;
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
            E_nonbond_ww.Ucoul += Vel;
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
            E_nonbond_ww.Ucoul += Vel;
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
            E_nonbond_ww.Ucoul += Vel;
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
            E_nonbond_ww.Ucoul += Vel;
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
            E_nonbond_ww.Ucoul += Vel;
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
            E_nonbond_ww.Ucoul += Vel;
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
            E_nonbond_ww.Ucoul += Vel;
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
            E_nonbond_ww.Ucoul += Vel;
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

    // printf("solvent: Ecoul = %f Evdw = %f\n", energies.Ucoul, energies.Uvdw);
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
            if (excluded[i] || excluded[j]) continue;
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

            E_nonbond_pw.Ucoul += Vela;
            E_nonbond_pw.Uvdw += (V_a - V_b);
        }
    }

    // printf("solute-solvent: Ecoul = %f Evdw = %f\n", energies.Ucoul, energies.Uvdw);
}