// TODO: Add impropers, bond pairs
#include "cpu/include/context.h"
#include "cpu/include/constants.h"

#include "system.h"
#include "qatoms.h"
#include "vdw_rules.h"
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
    auto &ctx = Context::instance();
    int i, j;
    coord_t da;
    double r2, r6, r;
    double ai_aii, aj_aii, ai_bii, aj_bii;
    q_catype_t qi_type;
    catype_t aj_type;
    bool bond23, bond14;
    double scaling, Vel, V_a, V_b, dv;

    for (int qi = 0; qi < ctx.n_qatoms; qi++) {
        for (int pj = 0; pj < ctx.n_patoms; pj++) {
            i = ctx.q_atoms[qi].a - 1;
            j = ctx.p_atoms[pj].a - 1;

            bond23 = ctx.LJ_matrix[i * ctx.n_atoms_solute + j] == 3;
            bond14 = ctx.LJ_matrix[i * ctx.n_atoms_solute + j] == 1;

            if (bond23) continue;
            if (ctx.excluded[i] || ctx.excluded[j]) continue;

            scaling = bond14 ? ctx.topo.el14_scale : 1;

            da.x = ctx.coords[j].x - ctx.coords[i].x;
            da.y = ctx.coords[j].y - ctx.coords[i].y;
            da.z = ctx.coords[j].z - ctx.coords[i].z;

            r2 = pow(da.x, 2) + pow(da.y, 2) + pow(da.z, 2);

            r6 = r2 * r2 * r2;
            r2 = 1 / r2;
            r = sqrt(r2);
            double r6inv = r2 * r2 * r2;  // 1/r^6 for vdW calculation

            for (int state = 0; state < ctx.n_lambdas; state++) {
                qi_type = ctx.q_catypes[ctx.q_atypes[qi + ctx.n_qatoms * state].code - 1];
                aj_type = ctx.catypes[ctx.atypes[j].code - 1];

                ai_aii = bond14 ? qi_type.Ai_14 : qi_type.Ai;
                aj_aii = bond14 ? aj_type.aii_1_4 : aj_type.aii_normal;
                ai_bii = bond14 ? qi_type.Bi_14 : qi_type.Bi;
                aj_bii = bond14 ? aj_type.bii_1_4 : aj_type.bii_normal;

                Vel = ctx.topo.coulomb_constant * scaling * ctx.q_charges[qi + ctx.n_qatoms * state].q * ctx.ccharges[ctx.charges[j].code - 1].charge * r;
                if (ctx.topo.vdw_rule == VDW_GEOMETRIC) {
                    calc_vdw_geometric(ai_aii, aj_aii, ai_bii, aj_bii, r6inv, &V_a, &V_b);
                } else {
                    calc_vdw_arithmetic(ai_aii, aj_aii, ai_bii, aj_bii, r6inv, &V_a, &V_b);
                }
                dv = r2 * (-Vel - (12 * V_a - 6 * V_b)) * ctx.lambdas[state];

                // Update forces
                ctx.dvelocities[i].x -= dv * da.x;
                ctx.dvelocities[i].y -= dv * da.y;
                ctx.dvelocities[i].z -= dv * da.z;
                ctx.dvelocities[j].x += dv * da.x;
                ctx.dvelocities[j].y += dv * da.y;
                ctx.dvelocities[j].z += dv * da.z;

                // Update Q totals
                ctx.EQ_nonbond_qp[state].Ucoul += Vel;
                ctx.EQ_nonbond_qp[state].Uvdw += (V_a - V_b);
            }
        }
    }
}


void calc_nonbonded_qw_forces() {
    auto &ctx = Context::instance();
    int i;
    coord_t dO, dH1, dH2;
    double r2O, rH1, rH2, r6O, rO, r2H1, r2H2;
    double dvO, dvH1, dvH2;
    double V_a, V_b, VelO, VelH1, VelH2;
    q_catype_t qi_type;
    double ai_aii, ai_bii;

    if (ctx.A_O == 0) {
        catype_t catype_ow;    // Atom type of first O, H atom

        catype_ow = ctx.catypes[ctx.atypes[ctx.n_atoms_solute].code - 1];

        ctx.A_O = catype_ow.aii_normal;
        ctx.B_O = catype_ow.bii_normal;
    }

    // Loop over O-atoms, q-atoms
    for (int j = ctx.n_atoms_solute; j < ctx.n_atoms; j+= 3) {
        for (int qi = 0; qi < ctx.n_qatoms; qi++) {
            i = ctx.q_atoms[qi].a - 1;
            if (ctx.excluded[i] || ctx.excluded[j]) continue;
            dO.x = ctx.coords[j].x - ctx.coords[i].x;
            dO.y = ctx.coords[j].y - ctx.coords[i].y;
            dO.z = ctx.coords[j].z - ctx.coords[i].z;
            dH1.x = ctx.coords[j+1].x - ctx.coords[i].x;
            dH1.y = ctx.coords[j+1].y - ctx.coords[i].y;
            dH1.z = ctx.coords[j+1].z - ctx.coords[i].z;
            dH2.x = ctx.coords[j+2].x - ctx.coords[i].x;
            dH2.y = ctx.coords[j+2].y - ctx.coords[i].y;
            dH2.z = ctx.coords[j+2].z - ctx.coords[i].z;
            r2O = pow(dO.x, 2) + pow(dO.y, 2) + pow(dO.z, 2);
            rH1 = sqrt(1.0 / (pow(dH1.x, 2) + pow(dH1.y, 2) + pow(dH1.z, 2)));
            rH2 = sqrt(1.0 / (pow(dH2.x, 2) + pow(dH2.y, 2) + pow(dH2.z, 2)));
            r6O = r2O * r2O * r2O;
            r2O = 1.0 / r2O;
            rO = sqrt(r2O);
            double r6Oinv = r2O * r2O * r2O;  // 1/r^6 for vdW calculation
            r2H1 = rH1 * rH1;
            r2H2 = rH2 * rH2;

            // Reset potential
            dvO = 0;
            dvH1 = 0;
            dvH2 = 0;

            for (int state = 0; state < ctx.n_lambdas; state++) {
                qi_type = ctx.q_catypes[ctx.q_atypes[qi + ctx.n_qatoms * state].code - 1];

                ai_aii = qi_type.Ai;
                ai_bii = qi_type.Bi;

                if (ctx.topo.vdw_rule == VDW_GEOMETRIC) {
                    calc_vdw_geometric(ai_aii, ctx.A_O, ai_bii, ctx.B_O, r6Oinv, &V_a, &V_b);
                } else {
                    calc_vdw_arithmetic(ai_aii, ctx.A_O, ai_bii, ctx.B_O, r6Oinv, &V_a, &V_b);
                }

                VelO = ctx.topo.coulomb_constant * ctx.crg_ow * ctx.q_charges[qi + ctx.n_qatoms * state].q * rO;
                VelH1 = ctx.topo.coulomb_constant * ctx.crg_hw * ctx.q_charges[qi + ctx.n_qatoms * state].q * rH1;
                VelH2 = ctx.topo.coulomb_constant * ctx.crg_hw * ctx.q_charges[qi + ctx.n_qatoms * state].q * rH2;

                // if (state == 0 && qi == 1) printf("j = %d ai__aii = %f A_O = %f B_O = %f V_a = %f V_b = %f r6O = %f\n", j, ai_aii, A_O, B_O, V_a, V_b, r6O);

                dvO += r2O * (-VelO - (12 * V_a - 6 * V_b)) * ctx.lambdas[state];
                dvH1 -= r2H1 * VelH1 * ctx.lambdas[state];
                dvH2 -= r2H2 * VelH2 * ctx.lambdas[state];

                ctx.EQ_nonbond_qw[state].Ucoul += (VelO + VelH1 + VelH2);
                ctx.EQ_nonbond_qw[state].Uvdw += (V_a - V_b);
            }

            // Note r6O is not the usual 1/rO^6, but rather rO^6. be careful!!!

            // Update forces on Q-atom
            ctx.dvelocities[i].x -= (dvO * dO.x + dvH1 * dH1.x + dvH2 * dH2.x);
            ctx.dvelocities[i].y -= (dvO * dO.y + dvH1 * dH1.y + dvH2 * dH2.y);
            ctx.dvelocities[i].z -= (dvO * dO.z + dvH1 * dH1.z + dvH2 * dH2.z);

            // Update forces on water
            ctx.dvelocities[j].x += dvO * dO.x;
            ctx.dvelocities[j].y += dvO * dO.y;
            ctx.dvelocities[j].z += dvO * dO.z;
            ctx.dvelocities[j+1].x += dvH1 * dH1.x;
            ctx.dvelocities[j+1].y += dvH1 * dH1.y;
            ctx.dvelocities[j+1].z += dvH1 * dH1.z;
            ctx.dvelocities[j+2].x += dvH2 * dH2.x;
            ctx.dvelocities[j+2].y += dvH2 * dH2.y;
            ctx.dvelocities[j+2].z += dvH2 * dH2.z;
        }
    }
    
    #ifdef DEBUG 
    printf("q-w: Ecoul = %f Evdw = %f\n", ctx.EQ_nonbond_qw[0].Ucoul, ctx.EQ_nonbond_qw[0].Uvdw);
    #endif
}


void calc_nonbonded_qq_forces() {
    auto &ctx = Context::instance();
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

    for (int state = 0; state < ctx.n_lambdas; state++) {
        for (int qi = 0; qi < ctx.n_qatoms; qi++) {
            for (int qj = qi+1; qj < ctx.n_qatoms; qj++) {
                ai = ctx.q_atoms[qi].a - 1;
                aj = ctx.q_atoms[qj].a - 1;

                crg_i = ctx.q_charges[qi + ctx.n_qatoms * state].q;
                crg_j = ctx.q_charges[qj + ctx.n_qatoms * state].q;

                bond23 = ctx.LJ_matrix[ai * ctx.n_atoms_solute + aj] == 3;
                bond14 = ctx.LJ_matrix[ai * ctx.n_atoms_solute + aj] == 1;
        
                if (bond23) continue;
                if (ctx.excluded[ai] || ctx.excluded[aj]) continue;
    
                scaling = bond14 ? ctx.topo.el14_scale : 1;

                elscale = 1;
                for (int k = 0; k < ctx.n_qelscales; k++) {
                    if (ctx.q_elscales[k + ctx.n_qelscales * state].qi == qi+1 && ctx.q_elscales[k + ctx.n_qelscales * state].qj == qj+1) {
                        elscale = ctx.q_elscales[k + ctx.n_qelscales * state].mu;
                    }
                }

                qi_type = ctx.q_catypes[ctx.q_atypes[qi + ctx.n_qatoms * state].code - 1];
                qj_type = ctx.q_catypes[ctx.q_atypes[qj + ctx.n_qatoms * state].code - 1];

                da.x = ctx.coords[aj].x - ctx.coords[ai].x;
                da.y = ctx.coords[aj].y - ctx.coords[ai].y;
                da.z = ctx.coords[aj].z - ctx.coords[ai].z;
                r2a = 1 / (pow(da.x, 2) + pow(da.y, 2) + pow(da.z, 2));
                ra = sqrt(r2a);
                r6a = r2a * r2a * r2a;

                Vela = scaling * ctx.topo.coulomb_constant * crg_i * crg_j * ra * elscale;

                ai_aii = bond14 ? qi_type.Ai_14 : qi_type.Ai;
                aj_aii = bond14 ? qj_type.Ai_14 : qj_type.Ai;
                ai_bii = bond14 ? qi_type.Bi_14 : qi_type.Bi;
                aj_bii = bond14 ? qj_type.Bi_14 : qj_type.Bi;

                if (ctx.topo.vdw_rule == VDW_GEOMETRIC) {
                    calc_vdw_geometric(ai_aii, aj_aii, ai_bii, aj_bii, r6a, &V_a, &V_b);
                } else {
                    calc_vdw_arithmetic(ai_aii, aj_aii, ai_bii, aj_bii, r6a, &V_a, &V_b);
                }
                dva = r2a * ( -Vela -12 * V_a + 6 * V_b) * ctx.lambdas[state];

                ctx.dvelocities[ai].x -= dva * da.x;
                ctx.dvelocities[ai].y -= dva * da.y;
                ctx.dvelocities[ai].z -= dva * da.z;

                ctx.dvelocities[aj].x += dva * da.x;
                ctx.dvelocities[aj].y += dva * da.y;
                ctx.dvelocities[aj].z += dva * da.z;
    
                ctx.EQ_nonbond_qq[state].Ucoul += Vela;
                ctx.EQ_nonbond_qq[state].Uvdw += (V_a - V_b);    
            }
        }
    }

    #ifdef DEBUG 
    printf("q-q: Ecoul = %f Evdw = %f\n", ctx.EQ_nonbond_qq[0].Ucoul, ctx.EQ_nonbond_qq[0].Uvdw);
    #endif
}

void calc_qangle_forces(int state) {
    auto &ctx = Context::instance();
    int ic;
    int ai, aj, ak;
    coord_t rji, rjk;
    double bji, bjk;
    double cos_th, th, dth, ener, dv, f1;
    coord_t di, dk;

    for (int i = 0; i < ctx.n_qangles; i++) {
        ic = ctx.q_angles[i + ctx.n_qangles * state].code-1;

        // Skip if angle not present (code 0)
        if (ic == 0) continue;

        ai = ctx.q_angles[i + ctx.n_qangles * state].ai - 1;
        aj = ctx.q_angles[i + ctx.n_qangles * state].aj - 1;
        ak = ctx.q_angles[i + ctx.n_qangles * state].ak - 1;

        rji.x = ctx.coords[ai].x - ctx.coords[aj].x;
        rji.y = ctx.coords[ai].y - ctx.coords[aj].y;
        rji.z = ctx.coords[ai].z - ctx.coords[aj].z;

        rjk.x = ctx.coords[ak].x - ctx.coords[aj].x;
        rjk.y = ctx.coords[ak].y - ctx.coords[aj].y;
        rjk.z = ctx.coords[ak].z - ctx.coords[aj].z;

        bji = sqrt(pow(rji.x, 2) + pow(rji.y, 2) + pow(rji.z, 2));
        bjk = sqrt(pow(rjk.x, 2) + pow(rjk.y, 2) + pow(rjk.z, 2));
        cos_th = rji.x * rjk.x + rji.y * rjk.y + rji.z * rjk.z;
        cos_th /= (bji * bjk);
        if (cos_th > 1) cos_th = 1;
        if (cos_th < -1) cos_th = -1;
        th = acos(cos_th);
        dth = th - to_radians(ctx.q_cangles[ic].th0);
        ener = .5 * ctx.q_cangles[ic].kth * pow(dth, 2);
        ctx.EQ_bond[state].Uangle += ener;

        dv = ctx.q_cangles[ic].kth * dth * ctx.lambdas[state];
        f1 = sin(th);
        if (abs(f1) < 1E-12) f1 = 1E-12;
        f1 = -1.0 / f1;

        di.x = f1 * (rjk.x / (bji * bjk) - cos_th * rji.x / pow(bji, 2));
        di.y = f1 * (rjk.y / (bji * bjk) - cos_th * rji.y / pow(bji, 2));
        di.z = f1 * (rjk.z / (bji * bjk) - cos_th * rji.z / pow(bji, 2));
        dk.x = f1 * (rji.x / (bji * bjk) - cos_th * rjk.x / pow(bjk, 2));
        dk.y = f1 * (rji.y / (bji * bjk) - cos_th * rjk.y / pow(bjk, 2));
        dk.z = f1 * (rji.z / (bji * bjk) - cos_th * rjk.z / pow(bjk, 2));

        ctx.dvelocities[ai].x += dv * di.x;
        ctx.dvelocities[ai].y += dv * di.y;
        ctx.dvelocities[ai].z += dv * di.z;
        ctx.dvelocities[ak].x += dv * dk.x;
        ctx.dvelocities[ak].y += dv * dk.y;
        ctx.dvelocities[ak].z += dv * dk.z;
        ctx.dvelocities[aj].x -= dv * (di.x + dk.x);
        ctx.dvelocities[aj].y -= dv * (di.y + dk.y);
        ctx.dvelocities[aj].z -= dv * (di.z + dk.z);
    }
}

void calc_qbond_forces(int state) {
    auto &ctx = Context::instance();
    int ic;
    int ai, aj;
    double b, db, ener, dv;
    coord_t rij;

    for (int i = 0; i < ctx.n_qbonds; i++) {
        ic = ctx.q_bonds[i + ctx.n_qbonds * state].code;

        if (ic == 0) continue;

        ai = ctx.q_bonds[i + ctx.n_qbonds * state].ai - 1;
        aj = ctx.q_bonds[i + ctx.n_qbonds * state].aj - 1;

        rij.x = ctx.coords[aj].x - ctx.coords[ai].x;
        rij.y = ctx.coords[aj].y - ctx.coords[ai].y;
        rij.z = ctx.coords[aj].z - ctx.coords[ai].z;

        b = sqrt(pow(rij.x, 2) + pow(rij.y, 2) + pow(rij.z, 2));
        db = b - ctx.q_cbonds[ic].b0;

        ener = 0.5 * ctx.q_cbonds[ic].kb * pow(db, 2);
        ctx.EQ_bond[state].Ubond += ener;
        dv = db * ctx.q_cbonds[ic].kb * ctx.lambdas[state] / b;

        ctx.dvelocities[ai].x -= dv * rij.x;
        ctx.dvelocities[ai].y -= dv * rij.y;
        ctx.dvelocities[ai].z -= dv * rij.z;
        ctx.dvelocities[aj].x += dv * rij.x;
        ctx.dvelocities[aj].y += dv * rij.y;
        ctx.dvelocities[aj].z += dv * rij.z;
    }
}

void calc_qtorsion_forces(int state) {
    auto &ctx = Context::instance();
    int ic;
    int ai, aj, ak, al;
    coord_t rji, rjk, rkl, rnj, rnk, rki, rlj;
    coord_t di, dl, dpi, dpj, dpk, dpl;

    double bj2inv, bk2inv, bjinv, bkinv;
    double bj, bk, cos_phi, phi;
    double arg, dv, f1;
    double ener;

    for (int i = 0; i < ctx.n_qtorsions; i++) {
        ic = ctx.q_torsions[i + ctx.n_qtorsions * state].code;

        if (ic == 0) continue;

        ai = ctx.q_torsions[i + ctx.n_qtorsions * state].ai - 1;
        aj = ctx.q_torsions[i + ctx.n_qtorsions * state].aj - 1;
        ak = ctx.q_torsions[i + ctx.n_qtorsions * state].ak - 1;
        al = ctx.q_torsions[i + ctx.n_qtorsions * state].al - 1;

        rji.x = ctx.coords[ai].x - ctx.coords[aj].x;
        rji.y = ctx.coords[ai].y - ctx.coords[aj].y;
        rji.z = ctx.coords[ai].z - ctx.coords[aj].z;
        rjk.x = ctx.coords[ak].x - ctx.coords[aj].x;
        rjk.y = ctx.coords[ak].y - ctx.coords[aj].y;
        rjk.z = ctx.coords[ak].z - ctx.coords[aj].z;
        rkl.x = ctx.coords[al].x - ctx.coords[ak].x;
        rkl.y = ctx.coords[al].y - ctx.coords[ak].y;
        rkl.z = ctx.coords[al].z - ctx.coords[ak].z;
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
        arg = ctx.q_ctorsions[ic].n * phi - to_radians(ctx.q_ctorsions[ic].d);
        ener = ctx.q_ctorsions[ic].k * (1 + cos(arg));
        dv = - ctx.q_ctorsions[ic].n * ctx.q_ctorsions[ic].k * sin(arg) * ctx.lambdas[state];

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
        ctx.EQ_bond[state].Utor += ener;

        ctx.dvelocities[ai].x += dv * dpi.x;
        ctx.dvelocities[ai].y += dv * dpi.y;
        ctx.dvelocities[ai].z += dv * dpi.z;

        ctx.dvelocities[aj].x += dv * dpj.x;
        ctx.dvelocities[aj].y += dv * dpj.y;
        ctx.dvelocities[aj].z += dv * dpj.z;

        ctx.dvelocities[ak].x += dv * dpk.x;
        ctx.dvelocities[ak].y += dv * dpk.y;
        ctx.dvelocities[ak].z += dv * dpk.z;

        ctx.dvelocities[al].x += dv * dpl.x;
        ctx.dvelocities[al].y += dv * dpl.y;
        ctx.dvelocities[al].z += dv * dpl.z;
    }
}
