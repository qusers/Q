#include "cpu_water_boundary_force.h"

#include "constants.h"
#include "cpu_force_accumulation.h"
#include "geometry.h"

void CpuWaterBoundaryForce::calc(Context& ctx, int iteration) {
    if (ctx.n_waters <= 0) return;
    calc_radix(ctx);
    if (ctx.md.polarisation) calc_polx(ctx, iteration);
}

void CpuWaterBoundaryForce::calc_radix(Context& ctx) {
    auto& coords = ctx.coords->cpu_data_p;
    auto& dvelocities = ctx.dvelocities->cpu_data_p;
    auto* excluded = ctx.excluded->cpu_data_p;
    const double* temperature_results = temperature_->data().results->cpu_data_p;  // was the param

    energy_accum_t uradx = 0;

    double shift;
    if (ctx.md.radial_force != 0.0) {
        shift = sqrt(Boltz * temperature_results[R_TFREE] / ctx.md.radial_force);
    } else {
        shift = 0.0;
    }

    for (int i = ctx.n_atoms_solute; i < ctx.n_atoms; i += 3) {
        if (excluded[i]) continue;
        const double dx = coords[i].x - ctx.topo.solvent_center.x;
        const double dy = coords[i].y - ctx.topo.solvent_center.y;
        const double dz = coords[i].z - ctx.topo.solvent_center.z;
        const double b = sqrt(dx * dx + dy * dy + dz * dz);
        const double db = b - (ctx.topo.solvent_radius - shift);

        double ener, dv;
        if (db > 0.0) {
            ener = 0.5 * ctx.md.radial_force * db * db - data_.Dwmz;  // ctx.Dwmz → data_.Dwmz
            dv = ctx.md.radial_force * db / b;
        } else if (b > 0.0) {
            const double fexp = exp(data_.awmz * db);  // ctx.awmz → data_.awmz
            ener = data_.Dwmz * (fexp * fexp - 2.0 * fexp);
            dv = -2.0 * data_.Dwmz * data_.awmz * (fexp - fexp * fexp) / b;
        } else {
            dv = 0.0;
            ener = 0.0;
        }
        add_energy(uradx, ener);
        add_force(dvelocities[i].x, dv * dx);
        add_force(dvelocities[i].y, dv * dy);
        add_force(dvelocities[i].z, dv * dz);
    }
    ctx.energy.host()[E_RESTR_RADX] += uradx;  // note: += (radix accumulates)
}

void CpuWaterBoundaryForce::calc_polx(Context& ctx, int iteration) {
    auto& coords = ctx.coords->cpu_data_p;
    auto& dvelocities = ctx.dvelocities->cpu_data_p;
    auto& exclude = ctx.excluded->cpu_data_p;

    auto* wshells = data_.wshells->cpu_data_p;

    auto* theta = data_.theta->cpu_data_p;
    auto* list_sh = data_.list_sh->cpu_data_p;

    const int n_shells = data_.n_shells;
    const int n_max_inshell = data_.n_max_inshell;

    if (iteration != 0 && iteration % itdis_update == 0) {
        for (int is = 0; is < n_shells; is++) {
            wshells[is].avtheta /= itdis_update;
            wshells[is].avn_inshell /= itdis_update;
            wshells[is].theta_corr = wshells[is].theta_corr + wshells[is].avtheta - acos(wshells[is].cstb);
            wshells[is].avtheta = 0.0;
            wshells[is].avn_inshell = 0.0;
        }
    }

    for (int is = 0; is < n_shells; is++) {
        wshells[is].n_inshell = 0;
    }

    for (int i = 0; i < ctx.n_waters; i++) {
        int wi = ctx.n_atoms_solute + 3 * i;
        if (exclude[wi]) continue;

        coord_t rmu = coords[wi + 1] + coords[wi + 2] - coords[wi] * 2;
        const double rm = norm(rmu);
        rmu = rmu / rm;

        coord_t rcu = coords[wi] - ctx.topo.solvent_center;
        const double rc = norm(rcu);
        rcu = rcu / rc;

        double cos_th = dot(rmu, rcu);
        cos_th = std::min(cos_th, 1.0);
        cos_th = std::max(cos_th, -1.0);
        theta[i] = acos(cos_th);

        if (rc > wshells[n_shells - 1].router - wshells[n_shells - 1].dr) {
            // In the shell
            int iis = n_shells - 1;
            while (iis > 0 && rc > wshells[iis].router) {
                iis--;
            }
            list_sh[iis * n_max_inshell + wshells[iis].n_inshell] = i;
            wshells[iis].n_inshell++;
        }
    }

    for (int i = 0; i < n_shells; i++) {
        int sz = wshells[i].n_inshell;
        if (sz == 0) continue;
        int offset = i * n_max_inshell;
        std::sort(list_sh + offset, list_sh + offset + sz, [&](int p, int q) {
            return theta[p] < theta[q];
        });
    }

    energy_accum_t uplox = 0;
    for (int is = 0; is < n_shells; is++) {
        int sz = wshells[is].n_inshell;
        if (sz == 0) continue;
        int offset = is * n_max_inshell;

        double avtdum = 0;
        for (int il = 0; il < sz; il++) {
            int ii = list_sh[offset + il];
            const double cos_arg = 1.0 + ((1.0 - 2.0 * (il + 1)) / sz);
            double arg = acos(cos_arg);
            arg = arg - 3.0 * sin(arg) * wshells[is].cstb / 2.0;
            arg = std::max(arg, 0.0);
            arg = std::min(arg, M_PI);

            avtdum += theta[ii];
            const double dtheta = theta[ii] - arg + wshells[is].theta_corr;
            const double ener = 0.5 * ctx.md.polarisation_force * dtheta * dtheta;

            add_energy(uplox, ener);

            const double dv = ctx.md.polarisation_force * dtheta;

            int wi = ctx.n_atoms_solute + 3 * ii;

            coord_t rmu = coords[wi + 1] + coords[wi + 2] - coords[wi] * 2;
            const double rm = norm(rmu);
            rmu = rmu / rm;

            coord_t rcu = coords[wi] - ctx.topo.solvent_center;
            const double rc = norm(rcu);
            rcu = rcu / rc;

            double cos_th = dot(rmu, rcu);
            cos_th = std::min(cos_th, 1.0);
            cos_th = std::max(cos_th, -1.0);

            double f0 = sin(acos(cos_th));
            if (std::abs(f0) < k_singular_sin_epsilon) f0 = k_singular_sin_epsilon;
            f0 = -dv / f0;

            // ∂cosθ/∂x, closed form (rmu, rcu already normalized; divide by original lengths rm, rc)
            const coord_t g_mu = (rcu - rmu * cos_th) / rm;  // dipole-vector side
            const coord_t g_c = (rmu - rcu * cos_th) / rc;   // radial-vector side

            const coord_t fO = g_mu * (-2.0) + g_c;  // O: ∂rmu/∂O=-2  and  ∂rcu/∂O=+1
            const coord_t fH = g_mu;                 // each H: ∂rmu/∂H=+1

            add_force(dvelocities[wi].x, f0 * fO.x);
            add_force(dvelocities[wi].y, f0 * fO.y);
            add_force(dvelocities[wi].z, f0 * fO.z);
            add_force(dvelocities[wi + 1].x, f0 * fH.x);
            add_force(dvelocities[wi + 1].y, f0 * fH.y);
            add_force(dvelocities[wi + 1].z, f0 * fH.z);
            add_force(dvelocities[wi + 2].x, f0 * fH.x);
            add_force(dvelocities[wi + 2].y, f0 * fH.y);
            add_force(dvelocities[wi + 2].z, f0 * fH.z);
        }
        wshells[is].avtheta += avtdum / wshells[is].n_inshell;
        wshells[is].avn_inshell += wshells[is].n_inshell;
    }
    ctx.energy.host()[E_RESTR_POLX] = uplox;
}