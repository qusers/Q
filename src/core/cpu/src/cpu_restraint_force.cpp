#include "cpu_restraint_force.h"

#include "bonded_force.h"
#include "constants.h"
#include "cpu_force_accumulation.h"
#include "geometry.h"
#include "q_math.h"

void CpuRestraintForce::calc(Context& ctx) {
    calc_posrestr(ctx);
    calc_restrseq(ctx);
    calc_restrpos(ctx);
    calc_restrdis(ctx);
    calc_restrang(ctx);
    calc_restrwall(ctx);
}

void CpuRestraintForce::calc_posrestr(Context& ctx) {
    if (data_.posrestr.n == 0) return;
    auto* ids = data_.posrestr.ids->cpu_data_p;
    auto* kk = data_.posrestr.k->cpu_data_p;
    auto* eslot_a = data_.posrestr.eslot_a->cpu_data_p;
    auto* eslot_b = data_.posrestr.eslot_b->cpu_data_p;
    auto* coords = ctx.coords->cpu_data_p;
    auto* coords_init = ctx.coords_init->cpu_data_p;
    auto* dvelocities = ctx.dvelocities->cpu_data_p;
    auto* e = ctx.energy.host();

    for (int i = 0; i < data_.posrestr.n; i++) {
        const int idx = ids[i];
        const double k = kk[i];

        const double dx = coords[idx].x - coords_init[idx].x;
        const double dy = coords[idx].y - coords_init[idx].y;
        const double dz = coords[idx].z - coords_init[idx].z;
        const double ener = 0.5 * k * (dx * dx + dy * dy + dz * dz);

        if (eslot_a[i] >= 0) add_energy(e[eslot_a[i]], ener);
        if (eslot_b[i] >= 0) add_energy(e[eslot_b[i]], ener);

        add_force(dvelocities[idx].x, k * dx);
        add_force(dvelocities[idx].y, k * dy);
        add_force(dvelocities[idx].z, k * dz);
    }
}

void CpuRestraintForce::calc_restrseq(Context& ctx) {
    if (data_.restrseq.n == 0) return;
    auto* recs = data_.restrseq.recs->cpu_data_p;
    auto* atypes = ctx.atypes->cpu_data_p;
    auto* catypes = ctx.catypes->cpu_data_p;
    auto* coords = ctx.coords->cpu_data_p;
    auto* coords_init = ctx.coords_init->cpu_data_p;
    auto* dvelocities = ctx.dvelocities->cpu_data_p;
    auto* heavy = ctx.heavy->cpu_data_p;

    energy_accum_t upres = 0;

    for (int s = 0; s < data_.restrseq.n; s++) {
        const double k = recs[s].k;
        coord_t d = {0, 0, 0};
        int n_ctr = 0;
        double totmass = 0;

        if (recs[s].to_center == 1) {
            for (int i = recs[s].ai - 1; i < recs[s].aj - 1; i++) {
                if (heavy[i] || recs[s].ih)  // heavy or include H
                {
                    n_ctr++;
                    d = d + coords[i] - coords_init[i];
                }
            }

            if (n_ctr > 0) {
                const double inv_n_ctr = 1.0 / static_cast<double>(n_ctr);
                d = d * inv_n_ctr;
                const double r2 = norm2(d);
                const double ener = 0.5 * k * r2;
                add_energy(upres, ener);

                for (int i = recs[s].ai - 1; i < recs[s].aj - 1; i++) {
                    if (heavy[i] || recs[s].ih) {
                        const double mass = catypes[atypes[i].code - 1].m;
                        const double tmp = mass / CARBON_MASS;  // how many carbon-masses is this atom. todo: why is it?

                        add_force(dvelocities[i].x, k * d.x * tmp);
                        add_force(dvelocities[i].y, k * d.y * tmp);
                        add_force(dvelocities[i].z, k * d.z * tmp);
                    }
                }
            }

        } else if (recs[s].to_center == 2) {
            for (int i = recs[s].ai - 1; i < recs[s].aj - 1; i++) {
                if (heavy[i] || recs[s].ih) {
                    const double mass = catypes[atypes[i].code - 1].m;
                    totmass += mass;
                    d = d + (coords[i] - coords_init[i]) * mass;
                }
            }

            if (totmass > 0) {
                d = d * (1.0 / totmass);
                const double r2 = norm2(d);

                const double ener = 0.5 * k * r2;
                add_energy(upres, ener);
                for (int i = recs[s].ai - 1; i < recs[s].aj - 1; i++) {
                    if (heavy[i] || recs[s].ih) {
                        add_force(dvelocities[i].x, k * d.x);
                        add_force(dvelocities[i].y, k * d.y);
                        add_force(dvelocities[i].z, k * d.z);
                    }
                }
            }
        } else {
            for (int i = recs[s].ai - 1; i < recs[s].aj - 1; i++) {
                if (heavy[i] || recs[s].ih) {
                    d = coords[i] - coords_init[i];
                    const double r2 = norm2(d);
                    const double ener = 0.5 * k * r2;
                    add_energy(upres, ener);

                    add_force(dvelocities[i].x, k * d.x);
                    add_force(dvelocities[i].y, k * d.y);
                    add_force(dvelocities[i].z, k * d.z);
                }
            }
        }
    }
    auto* energy = ctx.energy.host();
    energy[E_RESTR_PRES] += upres;
}

void CpuRestraintForce::calc_restrpos(Context& ctx) {
    if (data_.restrpos.n == 0) return;
    auto* recs = data_.restrpos.recs->cpu_data_p;
    auto* coords = ctx.coords->cpu_data_p;
    auto* dvelocities = ctx.dvelocities->cpu_data_p;
    auto* lambdas = ctx.lambdas->cpu_data_p;

    std::vector<energy_accum_t> urestr(ctx.n_lambdas(), 0);
    energy_accum_t upres = 0;

    for (int ir = 0; ir < data_.restrpos.n; ir++) {
        const int state = recs[ir].ipsi - 1;
        const int i = recs[ir].a - 1;

        coord_t d = coords[i] - recs[ir].x;

        const double lambda = (recs[ir].ipsi != 0) ? lambdas[state] : 1.0;

        coord_t k = recs[ir].k;
        const double ener = 0.5 * (k.x * d.x * d.x + k.y * d.y * d.y + k.z * d.z * d.z);

        add_force(dvelocities[i].x, k.x * d.x * lambda);
        add_force(dvelocities[i].y, k.y * d.y * lambda);
        add_force(dvelocities[i].z, k.z * d.z * lambda);

        if (recs[ir].ipsi == 0) {
            for (int k = 0; k < ctx.n_lambdas(); k++) {
                add_energy(urestr[k], ener);
            }

            if (ctx.n_lambdas() == 0) {
                add_energy(upres, ener);
            }
        } else {
            add_energy(urestr[state], ener);
        }
    }
    auto* energy = ctx.energy.host();
    for (int state = 0; state < ctx.n_lambdas(); state++) {
        energy[ctx.energy.eq_index(ENERGY_FIXED_COUNT, state, EQ_RESTR_URESTR)] += urestr[state];
    }
    energy[E_RESTR_PRES] += upres;
}

void CpuRestraintForce::calc_restrdis(Context& ctx) {
    if (data_.restrdis.n == 0) return;
    auto* recs = data_.restrdis.recs->cpu_data_p;
    auto* coords = ctx.coords->cpu_data_p;
    auto* dvelocities = ctx.dvelocities->cpu_data_p;
    auto* lambdas = ctx.lambdas->cpu_data_p;

    std::vector<energy_accum_t> urestr(ctx.n_lambdas(), 0);
    energy_accum_t upres = 0;

    for (int ir = 0; ir < data_.restrdis.n; ir++) {
        const int state = recs[ir].ipsi - 1;
        const int i = recs[ir].ai - 1;
        const int j = recs[ir].aj - 1;

        coord_t aij = coords[j] - coords[i];

        const double lambda = (recs[ir].ipsi != 0) ? lambdas[state] : 1.0;

        const double r = norm(aij);

        double dr;
        if (r < recs[ir].d1) {
            dr = r - recs[ir].d1;
        } else if (r > recs[ir].d2) {
            dr = r - recs[ir].d2;
        } else {
            continue;
        }

        const double ener = 0.5 * recs[ir].k * dr * dr;
        const double dv = lambda * recs[ir].k * dr / r;

        add_force(dvelocities[j].x, aij.x * dv);
        add_force(dvelocities[j].y, aij.y * dv);
        add_force(dvelocities[j].z, aij.z * dv);

        add_force(dvelocities[i].x, -aij.x * dv);
        add_force(dvelocities[i].y, -aij.y * dv);
        add_force(dvelocities[i].z, -aij.z * dv);

        if (recs[ir].ipsi == 0) {
            for (int k = 0; k < ctx.n_lambdas(); k++) {
                add_energy(urestr[k], ener);
            }

            if (ctx.n_lambdas() == 0) {
                add_energy(upres, ener);
            }

        } else {
            add_energy(urestr[state], ener);
        }
    }
    auto* energy = ctx.energy.host();
    for (int state = 0; state < ctx.n_lambdas(); state++) {
        energy[ctx.energy.eq_index(ENERGY_FIXED_COUNT, state, EQ_RESTR_URESTR)] += urestr[state];
    }
    energy[E_RESTR_PRES] += upres;
}

void CpuRestraintForce::calc_restrang(Context& ctx) {
    if (data_.restrang.n == 0) return;
    auto* recs = data_.restrang.recs->cpu_data_p;
    auto* coords = ctx.coords->cpu_data_p;
    auto* dvelocities = ctx.dvelocities->cpu_data_p;
    auto* lambdas = ctx.lambdas->cpu_data_p;

    std::vector<energy_accum_t> urestr(ctx.n_lambdas(), 0);
    energy_accum_t upres = 0;

    for (int ir = 0; ir < data_.restrang.n; ir++) {
        const int state = recs[ir].ipsi - 1;
        const double lambda = (recs[ir].ipsi != 0) ? lambdas[state] : 1.0;
        const int i = recs[ir].ai - 1;
        const int j = recs[ir].aj - 1;
        const int k = recs[ir].ak - 1;

        coord_t aji = coords[i] - coords[j];
        coord_t ajk = coords[k] - coords[j];

        double r_ji = norm(aji);
        double r_jk = norm(ajk);

        double cos_th = dot(aji, ajk) / (r_ji * r_jk);
        cos_th = std::min(cos_th, pt999);
        cos_th = std::max(cos_th, -pt999);
        const double th = acos(cos_th);

        // use calc_angle in bonded_force
        auto [v, dv_per_angle] = calc_angle(recs[ir].k, th, deg2rad(recs[ir].ang));
        coord_t perpendicular_v = get_perpendicular_vector(aji, ajk);
        double sin_th = sin(th);
        double f1 = sin_th;
        f1 = -1 / f1;
        perpendicular_v = perpendicular_v * f1 / r_ji;
        coord_t force1 = perpendicular_v * dv_per_angle * lambda;
        add_force(dvelocities[i].x, force1.x);
        add_force(dvelocities[i].y, force1.y);
        add_force(dvelocities[i].z, force1.z);

        coord_t perpendicular_v2 = get_perpendicular_vector(ajk, aji);
        perpendicular_v2 = perpendicular_v2 * f1 / r_jk;
        coord_t force2 = perpendicular_v2 * dv_per_angle * lambda;

        add_force(dvelocities[k].x, force2.x);
        add_force(dvelocities[k].y, force2.y);
        add_force(dvelocities[k].z, force2.z);

        add_force(dvelocities[j].x, -(force1.x + force2.x));
        add_force(dvelocities[j].y, -(force1.y + force2.y));
        add_force(dvelocities[j].z, -(force1.z + force2.z));

        if (recs[ir].ipsi == 0) {
            for (int lambda_idx = 0; lambda_idx < ctx.n_lambdas(); lambda_idx++) {
                add_energy(urestr[lambda_idx], v);
            }
            if (ctx.n_lambdas() == 0) {
                add_energy(upres, v);
            }
        } else {
            add_energy(urestr[state], v);
        }
    }
    auto* energy = ctx.energy.host();
    for (int state = 0; state < ctx.n_lambdas(); state++) {
        energy[ctx.energy.eq_index(ENERGY_FIXED_COUNT, state, EQ_RESTR_URESTR)] += urestr[state];
    }
    energy[E_RESTR_PRES] += upres;
}

void CpuRestraintForce::calc_restrwall(Context& ctx) {
    if (data_.restrwall.n == 0) return;
    auto* recs = data_.restrwall.recs->cpu_data_p;
    auto* coords = ctx.coords->cpu_data_p;
    auto* dvelocities = ctx.dvelocities->cpu_data_p;
    auto* heavy = ctx.heavy->cpu_data_p;

    energy_accum_t upres = 0;

    for (int ir = 0; ir < data_.restrwall.n; ir++) {
        const double k = recs[ir].k;
        for (int i = recs[ir].ai - 1; i < recs[ir].aj - 1; i++) {
            if (heavy[i] || recs[ir].ih) {
                coord_t d = coords[i] - ctx.topo.solvent_center;

                const double b = norm(d);
                const double db = b - recs[ir].d;
                double ener, dv;
                if (db > 0) {
                    ener = 0.5 * k * db * db - recs[ir].dMorse;
                    dv = k * db / b;
                } else {
                    const double fexp = exp(recs[ir].aMorse * db);
                    ener = recs[ir].dMorse * (fexp * fexp - 2.0 * fexp);
                    dv = -2.0 * recs[ir].dMorse * recs[ir].aMorse * (fexp - fexp * fexp) / b;
                }

                add_energy(upres, ener);

                add_force(dvelocities[i].x, dv * d.x);
                add_force(dvelocities[i].y, dv * d.y);
                add_force(dvelocities[i].z, dv * d.z);
            }
        }
    }
    auto* energy = ctx.energy.host();
    energy[EQ_RESTR_URESTR] += upres;
}
