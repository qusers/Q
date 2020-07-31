// TODO: Add impropers, bond pairs

#include "qatoms.h"
#include "system.h"
#include "utils.h"
#include <math.h>

void calc_nonbonded_qp_forces() {

}

void calc_nonbonded_qw_forces() {

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
        q_energies[state].Ucoul = 0;
        q_energies[state].Uvdw = 0;
        
        for (int qi = 0; qi < n_qatoms; qi++) {
            for (int qj = qi+1; qj < n_qatoms; qj++) {
                ai = q_atoms[qi].a - 1;
                aj = q_atoms[qj].a - 1;

                crg_i = q_charges[qi][state].q;
                crg_j = q_charges[qj][state].q;

                bond23 = false;
                bond14 = false;
                for (int k = 0; k < n_ngbrs23; k++) {
                    if (ngbrs23[k].ai == ai+1 && ngbrs23[k].aj == aj+1) {
                        // printf("BOND 23: %d, %d", ngbrs23[k].ai, ngbrs23[k].aj);
                        bond23 = true;
                    }
                }
                for (int k = 0; k < n_ngbrs14; k++) {
                    if (ngbrs14[k].ai == ai+1 && ngbrs14[k].aj == aj+1) {
                        // printf("BOND 14: %d, %d", ngbrs23[k].ai, ngbrs23[k].aj);
                        bond14 = true;
                    }
                }
    
                if (bond23) continue;
    
                scaling = bond14 ? .5 : 1;

                elscale = 0;
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