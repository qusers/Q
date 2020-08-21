// TODO: Add impropers

#include "system.h"
#include "bonded.h"
#include "utils.h"
#include <stdio.h>

/* =============================================
 * == BONDED INTERACTIONS
 * =============================================
 */

double calc_angle_forces(int start, int end) {
    int aii, aji, aki;

    coord_t ai, aj, ak;
    coord_t rji, rjk;
    coord_t di, dk;

    double bji2inv, bjk2inv, bjiinv, bjkinv;
    cangle_t cangle;
    double cos_th, th, dth, dv, f1;
    double ener;
    double angle = 0;

    for (int i = start; i < end; i++) {
        aii = angles[i].ai - 1;
        aji = angles[i].aj - 1;
        aki = angles[i].ak - 1;
        ai = coords[aii];
        aj = coords[aji];
        ak = coords[aki];

        cangle = cangles[angles[i].code - 1];

        rji.x = ai.x - aj.x;
        rji.y = ai.y - aj.y;
        rji.z = ai.z - aj.z;

        rjk.x = ak.x - aj.x;
        rjk.y = ak.y - aj.y;
        rjk.z = ak.z - aj.z;

        // Calculate inverse of norm of dist vector and their squares
        bji2inv = 1 / (rji.x * rji.x + rji.y * rji.y + rji.z * rji.z);
        bjk2inv = 1 / (rjk.x * rjk.x + rjk.y * rjk.y + rjk.z * rjk.z);
        bjiinv = sqrt(bji2inv);
        bjkinv = sqrt(bjk2inv);

        // Calculate cosine of angle and angle (th)
        cos_th = (rji.x * rjk.x + rji.y * rjk.y + rji.z * rjk.z) * bjiinv * bjkinv;

        if (cos_th > 1) {
            cos_th = 1;
        }
        else if (cos_th < -1) {
            cos_th = -1;
        }

        th = acos(cos_th);

        dth = th - to_radians(cangle.th0);
        ener = .5 * cangle.kth * pow(dth, 2);
        dv = cangle.kth*dth;

        f1 = sin(th);
        if (abs(f1) < 1.0E-12) {
            // Avoid division by zero
            f1 = -1.0E12;
        }
        else {
            f1 = -1.0 / f1;
        }

        // Update energies and forces

        angle += ener;

        di.x = f1 * ( rjk.x * bjiinv * bjkinv - cos_th * rji.x * bji2inv);
        di.y = f1 * ( rjk.y * bjiinv * bjkinv - cos_th * rji.y * bji2inv);
        di.z = f1 * ( rjk.z * bjiinv * bjkinv - cos_th * rji.z * bji2inv);

        dk.x = f1 * ( rji.x * bjiinv * bjkinv - cos_th * rjk.x * bjk2inv);
        dk.y = f1 * ( rji.y * bjiinv * bjkinv - cos_th * rjk.y * bjk2inv);
        dk.z = f1 * ( rji.z * bjiinv * bjkinv - cos_th * rjk.z * bjk2inv);

        dvelocities[aii].x += dv * di.x;
        dvelocities[aii].y += dv * di.y;
        dvelocities[aii].z += dv * di.z;

        dvelocities[aki].x += dv * dk.x;
        dvelocities[aki].y += dv * dk.y;
        dvelocities[aki].z += dv * dk.z;

        dvelocities[aji].x -= dv * (di.x + dk.x);
        dvelocities[aji].y -= dv * (di.y + dk.y);
        dvelocities[aji].z -= dv * (di.z + dk.z);

        // printf("ANGLE ener = %f\n", ener);
    }

    return angle;
}

double calc_bond_forces(int start, int end) {
    int aii, aji;
    coord_t ai, aj, dx;
    cbond_t cbond;
    double dx2, dx1, ddx, ener, ampl;
    double bond = 0;

    for (int i = start; i < end; i++) {
        aii = bonds[i].ai-1;
        aji = bonds[i].aj-1;
        ai = coords[aii];
        aj = coords[aji];

        cbond = cbonds[bonds[i].code-1];

        // Calculate distance vector, norm of distance vector
        dx.x = aj.x - ai.x;
        dx.y = aj.y - ai.y;
        dx.z = aj.z - ai.z;
        dx2 = dx.x * dx.x + dx.y * dx.y + dx.z * dx.z;
        dx1 = sqrt(dx2);

        // Calculate energy
        ddx = dx1 - cbond.b0;
        ener = .5 * cbond.kb * ddx * ddx;

        bond += ener;

        // Update forces
        ampl = cbond.kb * ddx / dx1;

        dvelocities[aji].x += ampl * dx.x;
        dvelocities[aji].y += ampl * dx.y;
        dvelocities[aji].z += ampl * dx.z;

        dvelocities[aii].x -= ampl * dx.x;
        dvelocities[aii].y -= ampl * dx.y;
        dvelocities[aii].z -= ampl * dx.z;

        // printf("BOND %d %d ener = %f\n", aii, aji, ener);
    }

    return bond;
}

double calc_torsion_forces(int start, int end) {
    int aii, aji, aki, ali;

    coord_t ai, aj, ak, al;
    coord_t rji, rjk, rkl, rnj, rnk, rki, rlj;
    coord_t di, dl, dpi, dpj, dpk, dpl;

    double bj2inv, bk2inv, bjinv, bkinv;
    double cos_phi, phi;
    double arg, dv, f1;
    double ener;
    double torsion = 0;

    torsion_t t;
    ctorsion_t ctors;

    for (int i = start; i < end; i++) {
        t = torsions[i];
        ctors = ctorsions[t.code - 1];

        aii = t.ai - 1;
        aji = t.aj - 1;
        aki = t.ak - 1;
        ali = t.al - 1;

        ai = coords[aii];
        aj = coords[aji];
        ak = coords[aki];
        al = coords[ali];

        rji.x = ai.x - aj.x;
        rji.y = ai.y - aj.y;
        rji.z = ai.z - aj.z;
        
        rjk.x = ak.x - aj.x;
        rjk.y = ak.y - aj.y;
        rjk.z = ak.z - aj.z;

        rkl.x = al.x - ak.x;
        rkl.y = al.y - ak.y;
        rkl.z = al.z - ak.z;

        rnj.x = rji.y * rjk.z - rji.z * rjk.y;
        rnj.y = rji.z * rjk.x - rji.x * rjk.z;
        rnj.z = rji.x * rjk.y - rji.y * rjk.x;

        rnk.x = -rjk.y * rkl.z + rjk.z * rkl.y;
        rnk.y = -rjk.z * rkl.x + rjk.x * rkl.z;
        rnk.z = -rjk.x * rkl.y + rjk.y * rkl.x;

        bj2inv = 1 / (pow(rnj.x, 2) + pow(rnj.y, 2) + pow(rnj.z, 2));
        bk2inv = 1 / (pow(rnk.x, 2) + pow(rnk.y, 2) + pow(rnk.z, 2));
        bjinv = sqrt(bj2inv);
        bkinv = sqrt(bk2inv);

        cos_phi = (rnj.x * rnk.x + rnj.y * rnk.y + rnj.z * rnk.z) * (bjinv * bkinv);
        if (cos_phi > 1) {
            cos_phi = 1;
            phi = acos(1);
        }
        else if (cos_phi < -1) {
            cos_phi = -1;
            phi = acos(-1);
        }
        else {
            phi = acos(cos_phi);
        }
        if (rjk.x * (rnj.y * rnk.z - rnj.z * rnk.y) 
                + rjk.y * (rnj.z * rnk.x - rnj.x * rnk.z)
                + rjk.z * (rnj.x * rnk.y - rnj.y * rnk.x) < 0) {
            phi = -phi;
        }

        
        // Energy
        arg = ctors.n * phi - to_radians(ctors.d);
        ener = ctors.k * (1 + cos(arg));
        dv = - ctors.n * ctors.k * sin(arg);

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
        torsion += ener;

        dvelocities[aii].x += dv * dpi.x;
        dvelocities[aii].y += dv * dpi.y;
        dvelocities[aii].z += dv * dpi.z;

        dvelocities[aji].x += dv * dpj.x;
        dvelocities[aji].y += dv * dpj.y;
        dvelocities[aji].z += dv * dpj.z;

        dvelocities[aki].x += dv * dpk.x;
        dvelocities[aki].y += dv * dpk.y;
        dvelocities[aki].z += dv * dpk.z;

        dvelocities[ali].x += dv * dpl.x;
        dvelocities[ali].y += dv * dpl.y;
        dvelocities[ali].z += dv * dpl.z;
    }

    return torsion;
}

double calc_improper2_forces(int start, int end) {
    int aii, aji, aki, ali;

    coord_t ai, aj, ak, al;
    coord_t rji, rjk, rkl, rnj, rnk, rki, rlj;
    double bj2inv, bk2inv, bjinv, bkinv;
    double cos_phi, phi, arg, ener, dv, f1;
    coord_t di, dl, dpi, dpj, dpk, dpl; 

    improper_t imp;
    cimproper_t cimp;
    double improper = 0;

    for (int i = start; i < end; i++) {
        imp = impropers[i];
        cimp = cimpropers[imp.code - 1];

        aii = imp.ai - 1;
        aji = imp.aj - 1;
        aki = imp.ak - 1;
        ali = imp.al - 1;

        ai = coords[aii];
        aj = coords[aji];
        ak = coords[aki];
        al = coords[ali];

        rji.x = ai.x - aj.x;
        rji.y = ai.y - aj.y;
        rji.z = ai.z - aj.z;
        rjk.x = ak.x - aj.x;
        rjk.y = ak.y - aj.y;
        rjk.z = ak.z - aj.z;
        rkl.x = al.x - ak.x;
        rkl.y = al.y - ak.y;
        rkl.z = al.z - ak.z;
        rnj.x = rji.y * rjk.z - rji.z * rjk.y;
        rnj.y = rji.z * rjk.x - rji.x * rjk.z;
        rnj.z = rji.x * rjk.y - rji.y * rjk.x;
        rnk.x = -rjk.y * rkl.z + rjk.z * rkl.y;
        rnk.y = -rjk.z * rkl.x + rjk.x * rkl.z;
        rnk.z = -rjk.x * rkl.y + rjk.y * rkl.x;

        bj2inv = 1 / (pow(rnj.x, 2) + pow(rnj.y, 2) + pow(rnj.z, 2));
        bk2inv = 1 / (pow(rnk.x, 2) + pow(rnk.y, 2) + pow(rnk.z, 2));
        bjinv = sqrt(bj2inv);
        bkinv = sqrt(bk2inv);

        cos_phi = (rnj.x * rnk.x + rnj.y * rnk.y + rnj.z * rnk.z) * (bjinv * bkinv);
        // printf("cos_phi = %f\n", cos_phi);
        if (cos_phi > 1) {
            cos_phi = 1;
        }
        if (cos_phi < -1) {
            cos_phi = -1;
        }
        phi = acos(cos_phi);
        if (rjk.x * (rnj.y * rnk.z - rnj.z * rnk.y)
            + rjk.y * (rnj.z * rnk.x - rnj.x * rnk.z)
            + rjk.z * (rnj.x * rnk.y - rnj.y * rnk.x) < 0) {
            phi = -phi;
        }

        // Energy
        arg = 2 * phi - to_radians(cimp.phi0);
        ener = cimp.k * (1 + cos(arg));
        dv = -2 * cimp.k * sin(arg);

        // Forces
        f1 = sin(phi);
        if (abs(f1) < 1E-12) f1 = 1E-12;
        f1 = -1 / f1;
        // printf("f1 = %f phi = %f cos_phi = %f\n", f1, phi, cos_phi);

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
        improper += ener;

        dvelocities[aii].x += dv * dpi.x;
        dvelocities[aii].y += dv * dpi.y;
        dvelocities[aii].z += dv * dpi.z;

        dvelocities[aji].x += dv * dpj.x;
        dvelocities[aji].y += dv * dpj.y;
        dvelocities[aji].z += dv * dpj.z;

        dvelocities[aki].x += dv * dpk.x;
        dvelocities[aki].y += dv * dpk.y;
        dvelocities[aki].z += dv * dpk.z;

        dvelocities[ali].x += dv * dpl.x;
        dvelocities[ali].y += dv * dpl.y;
        dvelocities[ali].z += dv * dpl.z;
    }

    return improper;
}