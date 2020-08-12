// TODO: Add impropers, bond pairs

#include "qatoms.h"
#include "system.h"
#include "utils.h"
#include <math.h>
#include <stdio.h>

void calc_nonbonded_qp_forces() {
    int i, j;
    coord_t da;
    double r2, r6, r;
    double ai_aii, aj_aii, ai_bii, aj_bii;
    q_catype_t qi_type;
    catype_t aj_type;
    bool bond23, bond14;
    double scaling, Vel, V_a, V_b, dv;

    // printf("n_qatoms = %d n_patoms = %d\n", n_qatoms, n_patoms);

    for (int qi = 0; qi < n_qatoms; qi++) {
        for (int pj = 0; pj < n_patoms; pj++) {
            i = q_atoms[qi].a - 1;
            j = p_atoms[pj].a - 1;

            printf("i = %d j = %d\n", i, j);

            bond23 = LJ_matrix[i * n_atoms_solute + j] == 3;
            bond14 = LJ_matrix[i * n_atoms_solute + j] == 1;

            if (bond23) continue;
            if (excluded[i] || excluded[j]) continue;

            scaling = bond14 ? .5 : 1;

            da.x = coords[j].x - coords[i].x;
            da.y = coords[j].y - coords[i].y;
            da.z = coords[j].z - coords[i].z;

            r2 = pow(da.x, 2) + pow(da.y, 2) + pow(da.z, 2);
            
            r6 = r2 * r2 * r2;
            r2 = 1 / r2;
            r = sqrt(r2);

            for (int state = 0; state < n_lambdas; state++) {
                qi_type = q_catypes[q_atypes[qi][state].code - 1];
                aj_type = catypes[atypes[j].code - 1];
    
                ai_aii = bond14 ? qi_type.Ai_14 : qi_type.Ai;
                aj_aii = bond14 ? aj_type.aii_1_4 : aj_type.aii_normal;
                ai_bii = bond14 ? qi_type.Bi_14 : qi_type.Bi;
                aj_bii = bond14 ? aj_type.bii_1_4 : aj_type.bii_normal;
    
                Vel = Coul * scaling * q_charges[qi][state].q * ccharges[charges[j].code - 1].charge * r;
                V_a = ai_aii * aj_aii / (r6 * r6);
                V_b = ai_bii * aj_bii / r6;
                dv = r2 * (-Vel - (12 * V_a - 6 * V_b)) * lambdas[state];

                // printf("Vel = %f V_a = %f V_b = %f\n", Vel, V_a, V_b);
    
                // Update forces
                dvelocities[i].x -= dv * da.x;
                dvelocities[i].y -= dv * da.y;
                dvelocities[i].z -= dv * da.z;
                dvelocities[j].x += dv * da.x;
                dvelocities[j].y += dv * da.y;
                dvelocities[j].z += dv * da.z;
    
                // Update Q totals
                q_energies[state].Ucoul += Vel;
                q_energies[state].Uvdw += (V_a - V_b);
            }
        }
    }

    // printf("q_energies[0].Ucoul = %f\n", q_energies[0].Ucoul);
    // printf("q_energies[0].Uvdw = %f\n", q_energies[0].Uvdw);
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
                qi_type = q_catypes[q_atypes[qi][state].code - 1];

                ai_aii = qi_type.Ai;
                ai_bii = qi_type.Bi;
                
                V_a = ai_aii * A_O / (r6O * r6O);
                V_b = ai_bii * B_O / (r6O);

                VelO = Coul * crg_ow * q_charges[qi][state].q * rO;
                VelH1 = Coul * crg_hw * q_charges[qi][state].q * rH1;
                VelH2 = Coul * crg_hw * q_charges[qi][state].q * rH2;
                dvO += r2O * (-VelO - (12 * V_a - 6 * V_b)) * lambdas[state];
                dvH1 -= r2H1 * VelH1 * lambdas[state];
                dvH2 -= r2H2 * VelH2 * lambdas[state];

                q_energies[state].Ucoul += (VelO + VelH1 + VelH2);
                q_energies[state].Uvdw += (V_a - V_b);
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

                crg_i = q_charges[qi][state].q;
                crg_j = q_charges[qj][state].q;

                bond23 = LJ_matrix[ai * n_atoms_solute + aj] == 3;
                bond14 = LJ_matrix[ai * n_atoms_solute + aj] == 1;
        
                if (bond23) continue;
                if (excluded[ai] || excluded[aj]) continue;
    
                scaling = bond14 ? .5 : 1;

                elscale = 1;
                for (int k = 0; k < n_qelscales; k++) {
                    if (q_elscales[k][state].qi == qi+1 && q_elscales[k][state].qj == qj+1) {
                        elscale = q_elscales[k][state].mu;
                    }
                }

                qi_type = q_catypes[q_atypes[qi][state].code - 1];
                qj_type = q_catypes[q_atypes[qj][state].code - 1];

                da.x = coords[aj].x - coords[ai].x;
                da.y = coords[aj].y - coords[ai].y;
                da.z = coords[aj].z - coords[ai].z;
                r2a = 1 / (pow(da.x, 2) + pow(da.y, 2) + pow(da.z, 2));
                ra = sqrt(r2a);
                r6a = r2a * r2a * r2a;
    
                Vela = scaling * Coul * crg_i * crg_j * ra * elscale;
    
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
    
                q_energies[state].Ucoul += Vela;
                q_energies[state].Uvdw += (V_a - V_b);    
            }
        }
    }
}

void calc_qangle_forces(int state) {
    int ic;
    int ai, aj, ak;
    coord_t rji, rjk;
    double bji, bjk;
    double cos_th, th, dth, ener, dv, f1;
    coord_t di, dk;

    for (int i = 0; i < n_qangles; i++) {
        ic = q_angles[i][state].code-1;

        // Skip if angle not present (code 0)
        if (ic == 0) continue;

        ai = q_angles[i][state].ai - 1;
        aj = q_angles[i][state].aj - 1;
        ak = q_angles[i][state].ak - 1;

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
        q_energies[state].Uangle += ener;

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
        ic = q_bonds[i][state].code;

        if (ic == 0) continue;

        ai = q_bonds[i][state].ai - 1;
        aj = q_bonds[i][state].aj - 1;

        rij.x = coords[aj].x - coords[ai].x;
        rij.y = coords[aj].y - coords[ai].y;
        rij.z = coords[aj].z - coords[ai].z;

        b = sqrt(pow(rij.x, 2) + pow(rij.y, 2) + pow(rij.z, 2));
        db = b - q_cbonds[ic].b0;

        ener = 0.5 * q_cbonds[ic].kb * pow(db, 2);
        q_energies[state].Ubond += ener;
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
        ic = q_torsions[i][state].code;

        if (ic == 0) continue;

        ai = q_torsions[i][state].ai - 1;
        aj = q_torsions[i][state].aj - 1;
        ak = q_torsions[i][state].ak - 1;
        al = q_torsions[i][state].al - 1;

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
        q_energies[state].Utor += ener;

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