#include "cpu_restrang_force.h"

#include <math.h>
#include <vector>

#include "context.h"
#include "cpu_force_accumulation.h"
#include "cpu_utils.h"

void calc_restrang_forces() {
    auto& ctx = Context::instance();
    auto &coords = ctx.coords->cpu_data_p;
    auto &dvelocities = ctx.dvelocities->cpu_data_p;
    auto &restrangs = ctx.restrangs->cpu_data_p;
    auto *lambdas = ctx.lambdas->cpu_data_p;
    auto *EQ_restraint = ctx.EQ_restraint->cpu_data_p;

    int state, i, j, k;
    std::vector<energy_accum_t> urestr(ctx.n_lambdas, 0);
    energy_accum_t upres = 0;

    for (int ir = 0; ir < ctx.n_restrangs; ir++) {
        state = restrangs[ir].ipsi - 1;
        i = restrangs[ir].ai - 1;
        j = restrangs[ir].aj - 1;
        k = restrangs[ir].ak - 1;

        const double dr_x = static_cast<double>(coords[i].x) - static_cast<double>(coords[j].x);
        const double dr_y = static_cast<double>(coords[i].y) - static_cast<double>(coords[j].y);
        const double dr_z = static_cast<double>(coords[i].z) - static_cast<double>(coords[j].z);

        const double dr2_x = static_cast<double>(coords[k].x) - static_cast<double>(coords[j].x);
        const double dr2_y = static_cast<double>(coords[k].y) - static_cast<double>(coords[j].y);
        const double dr2_z = static_cast<double>(coords[k].z) - static_cast<double>(coords[j].z);

        double lambda;
        if (restrangs[ir].ipsi != 0) {
            lambda = static_cast<double>(lambdas[state]);
        } else {
            lambda = 1.0;
        }

        const double r2ij = dr_x * dr_x + dr_y * dr_y + dr_z * dr_z;
        const double r2jk = dr2_x * dr2_x + dr2_y * dr2_y + dr2_z * dr2_z;

        const double rij = sqrt(r2ij);
        const double rjk = sqrt(r2jk);

        double cos_th = dr_x * dr2_x + dr_y * dr2_y + dr_z * dr2_z;
        cos_th /= rij * rjk;

        if (cos_th > 1) {
            cos_th = 1;
        }
        if (cos_th < -1) {
            cos_th = -1;
        }

        const double th = acos(cos_th);
        const double dth = th - to_radians(restrangs[ir].ang);

        const double ener = 0.5 * restrangs[ir].k * dth * dth;
        const double dv = lambda * restrangs[ir].k * dth;

        double f1 = sin(th);
        const double sin_epsilon = static_cast<double>(k_singular_sin_epsilon);
        if (fabs(f1) < sin_epsilon) {
            f1 = -1.0 / sin_epsilon;
        } else {
            f1 = -1 / f1;
        }

        const double di_x = f1 * (dr2_x / (rij * rjk) - cos_th * dr_x / r2ij);
        const double di_y = f1 * (dr2_y / (rij * rjk) - cos_th * dr_y / r2ij);
        const double di_z = f1 * (dr2_z / (rij * rjk) - cos_th * dr_z / r2ij);
        const double dk_x = f1 * (dr_x / (rij * rjk) - cos_th * dr2_x / r2jk);
        const double dk_y = f1 * (dr_y / (rij * rjk) - cos_th * dr2_y / r2jk);
        const double dk_z = f1 * (dr_z / (rij * rjk) - cos_th * dr2_z / r2jk);

        add_force(dvelocities[i].x, dv * di_x);
        add_force(dvelocities[i].y, dv * di_y);
        add_force(dvelocities[i].z, dv * di_z);
        add_force(dvelocities[k].x, dv * dk_x);
        add_force(dvelocities[k].y, dv * dk_y);
        add_force(dvelocities[k].z, dv * dk_z);
        add_force(dvelocities[j].x, -dv * (di_x + dk_x));
        add_force(dvelocities[j].y, -dv * (di_y + dk_y));
        add_force(dvelocities[j].z, -dv * (di_z + dk_z));

        if (restrangs[ir].ipsi == 0) {
            for (int lambda_idx = 0; lambda_idx < ctx.n_lambdas; lambda_idx++) {
                add_energy(urestr[lambda_idx], ener);
            }
            if (ctx.n_lambdas == 0) {
                add_energy(upres, ener);
            }
        } else {
            add_energy(urestr[state], ener);
        }
    }

    for (int state = 0; state < ctx.n_lambdas; state++) {
        EQ_restraint[state].Urestr += energy_from_accum(urestr[state]);
    }
    ctx.E_restraint.Upres += energy_from_accum(upres);
}
