#include "cpu_restrang_force.h"

#include <math.h>

#include "context.h"
#include "cpu_utils.h"

void calc_restrang_forces() {
    auto& ctx = Context::instance();

    int state, i, j, k;
    coord_t dr, dr2, di, dk;
    double lambda, r2ij, r2jk, rij, rjk, cos_th, th;
    double dth, dv, ener, f1;

    for (int ir = 0; ir < ctx.n_restrangs; ir++) {
        state = ctx.restrangs[ir].ipsi - 1;
        i = ctx.restrangs[ir].ai - 1;
        j = ctx.restrangs[ir].aj - 1;
        k = ctx.restrangs[ir].ak - 1;

        dr.x = ctx.coords[i].x - ctx.coords[j].x;
        dr.y = ctx.coords[i].y - ctx.coords[j].y;
        dr.z = ctx.coords[i].z - ctx.coords[j].z;

        dr2.x = ctx.coords[k].x - ctx.coords[j].x;
        dr2.y = ctx.coords[k].y - ctx.coords[j].y;
        dr2.z = ctx.coords[k].z - ctx.coords[j].z;

        if (ctx.restrangs[ir].ipsi != 0) {
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
        dth = th - to_radians(ctx.restrangs[ir].ang);

        ener = .5 * ctx.restrangs[ir].k * pow(dth, 2);
        dv = lambda * ctx.restrangs[ir].k * dth;

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

        ctx.dvelocities[i].x += dv * di.x;
        ctx.dvelocities[i].y += dv * di.y;
        ctx.dvelocities[i].z += dv * di.z;
        ctx.dvelocities[k].x += dv * dk.x;
        ctx.dvelocities[k].y += dv * dk.y;
        ctx.dvelocities[k].z += dv * dk.z;
        ctx.dvelocities[j].x -= dv * (di.x + dk.x);
        ctx.dvelocities[j].y -= dv * (di.y + dk.y);
        ctx.dvelocities[j].z -= dv * (di.z + dk.z);

        if (ctx.restrangs[ir].ipsi == 0) {
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
