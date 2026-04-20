#include "cpu_restrang_force.h"

#include <math.h>

#include "context.h"
#include "cpu_utils.h"

void calc_restrang_forces() {
    auto& ctx = Context::instance();
    auto &coords = ctx.coords->cpu_data_p;
    auto &dvelocities = ctx.dvelocities->cpu_data_p;
    auto &restrangs = ctx.restrangs->cpu_data_p;

    int state, i, j, k;
    coord_t dr, dr2, di, dk;
    double lambda, r2ij, r2jk, rij, rjk, cos_th, th;
    double dth, dv, ener, f1;

    for (int ir = 0; ir < ctx.n_restrangs; ir++) {
        state = restrangs[ir].ipsi - 1;
        i = restrangs[ir].ai - 1;
        j = restrangs[ir].aj - 1;
        k = restrangs[ir].ak - 1;

        dr.x = coords[i].x - coords[j].x;
        dr.y = coords[i].y - coords[j].y;
        dr.z = coords[i].z - coords[j].z;

        dr2.x = coords[k].x - coords[j].x;
        dr2.y = coords[k].y - coords[j].y;
        dr2.z = coords[k].z - coords[j].z;

        if (restrangs[ir].ipsi != 0) {
            lambda = ctx.lambdas[state];
        } else {
            lambda = 1;
        }

        r2ij = pow(dr.x, 2) + pow(dr.y, 2) + pow(dr.z, 2);
        r2jk = pow(dr2.x, 2) + pow(dr2.y, 2) + pow(dr2.z, 2);

        rij = sqrt(r2ij);
        rjk = sqrt(r2jk);

        cos_th = dr.x * dr2.x + dr.y * dr2.y + dr.z * dr2.z;
        cos_th /= rij * rjk;

        if (cos_th > 1) {
            cos_th = 1;
        }
        if (cos_th < -1) {
            cos_th = -1;
        }

        th = acos(cos_th);
        dth = th - to_radians(restrangs[ir].ang);

        ener = .5 * restrangs[ir].k * pow(dth, 2);
        dv = lambda * restrangs[ir].k * dth;

        f1 = sin(th);
        if (fabs(f1) < 1E-12) {
            f1 = -1E-12;
        } else {
            f1 = -1 / f1;
        }

        di.x = f1 * (dr2.x / (rij * rjk) - cos_th * dr.x / r2ij);
        di.y = f1 * (dr2.y / (rij * rjk) - cos_th * dr.y / r2ij);
        di.z = f1 * (dr2.z / (rij * rjk) - cos_th * dr.z / r2ij);
        dk.x = f1 * (dr.x / (rij * rjk) - cos_th * dr2.x / r2jk);
        dk.y = f1 * (dr.y / (rij * rjk) - cos_th * dr2.y / r2jk);
        dk.z = f1 * (dr.z / (rij * rjk) - cos_th * dr2.z / r2jk);

        dvelocities[i].x += dv * di.x;
        dvelocities[i].y += dv * di.y;
        dvelocities[i].z += dv * di.z;
        dvelocities[k].x += dv * dk.x;
        dvelocities[k].y += dv * dk.y;
        dvelocities[k].z += dv * dk.z;
        dvelocities[j].x -= dv * (di.x + dk.x);
        dvelocities[j].y -= dv * (di.y + dk.y);
        dvelocities[j].z -= dv * (di.z + dk.z);

        if (restrangs[ir].ipsi == 0) {
            for (int lambda_idx = 0; lambda_idx < ctx.n_lambdas; lambda_idx++) {
                ctx.EQ_restraint[lambda_idx].Urestr += ener;
            }
            if (ctx.n_lambdas == 0) {
                ctx.E_restraint.Upres += ener;
            }
        } else {
            ctx.EQ_restraint[state].Urestr += ener;
        }
    }
}
