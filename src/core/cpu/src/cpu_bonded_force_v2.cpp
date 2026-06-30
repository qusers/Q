#include "cpu_bonded_force_v2.h"

#include "constants.h"
#include "cpu_force_accumulation.h"
#include "geometry.h"

void CpuBondedForce::calc(Context& ctx) {
    calc_bonds(ctx);
    calc_angles(ctx);
    calc_torsions(ctx);
    calc_impropers(ctx);
}

void CpuBondedForce::calc_bonds(Context& ctx) {
    auto* ids = data_.bond.ids->cpu_data_p;
    auto* params = data_.bond.params->cpu_data_p;
    auto* eslot = data_.bond.eslot->cpu_data_p;
    auto* dvelocities = ctx.dvelocities->cpu_data_p;
    auto* coords = ctx.coords->cpu_data_p;

    energy_accum_t e[2] = {0, 0};  // [bonded_p_slot(), bonded_w_slot()]

    for (int i = 0; i < data_.bond.n; i++) {
        const int ai = ids[i].i, aj = ids[i].j;
        const double k = params[i].x, r_eq = params[i].y;

        const double dx = coords[aj].x - coords[ai].x;  // i -> j
        const double dy = coords[aj].y - coords[ai].y;
        const double dz = coords[aj].z - coords[ai].z;
        const double r = sqrt(dx * dx + dy * dy + dz * dz);

        auto [v, dv] = calc_bond(k, r, r_eq);
        // force = dv * unit_vector, where unit_vector = (dx,dy,dz)/r.
        // Fold the 1/r into ampl so we can multiply by the raw (dx,dy,dz).

        const double ampl = dv / r;  //
        const double fx = ampl * dx;
        const double fy = ampl * dy;
        const double fz = ampl * dz;

        add_force(dvelocities[ai].x, -fx);
        add_force(dvelocities[ai].y, -fy);
        add_force(dvelocities[ai].z, -fz);

        add_force(dvelocities[aj].x, fx);
        add_force(dvelocities[aj].y, fy);
        add_force(dvelocities[aj].z, fz);

        add_energy(e[eslot[i]], v);
    }

    ctx.E_bond_p.Ubond = energy_from_accum(e[bonded_p_slot()]);
    ctx.E_bond_w.Ubond = energy_from_accum(e[bonded_w_slot()]);
}

void CpuBondedForce::calc_angles(Context& ctx) {
    auto* ids = data_.angle.ids->cpu_data_p;
    auto* params = data_.angle.params->cpu_data_p;
    auto* eslot = data_.angle.eslot->cpu_data_p;
    auto* dvelocities = ctx.dvelocities->cpu_data_p;
    auto* coords = ctx.coords->cpu_data_p;

    energy_accum_t e[2] = {0, 0};  // [bonded_p_slot(), bonded_w_slot()]

    for (int i = 0; i < data_.angle.n; i++) {
        const int ai = ids[i].i, aj = ids[i].j, ak = ids[i].k;
        const double k = params[i].x;
        const double th_eq = params[i].y;

        coord_t aji = coords[ai] - coords[aj];
        coord_t ajk = coords[ak] - coords[aj];

        double r_ji = norm(aji);
        double r_jk = norm(ajk);

        double cos_th = dot(aji, ajk) / (r_ji * r_jk);
        cos_th = std::min(cos_th, pt999);
        cos_th = std::max(cos_th, -pt999);
        const double th = acos(cos_th);
        auto [v, dv_per_angle] = calc_angle(k, th, th_eq);

        add_energy(e[eslot[i]], v);

        // Now we should get a vector that is perpendicular to aji inside the plane formed by atoms i-j-k
        coord_t perpendicular_v = get_perpendicular_vector(aji, ajk);
        double sin_th = sin(th);
        double f1 = sin_th;
        f1 = -1.0 / f1;
        perpendicular_v = perpendicular_v * f1 / r_ji;
        coord_t force1 = perpendicular_v * dv_per_angle;
        add_force(dvelocities[ai].x, force1.x);
        add_force(dvelocities[ai].y, force1.y);
        add_force(dvelocities[ai].z, force1.z);

        coord_t perpendicular_v2 = get_perpendicular_vector(ajk, aji);
        perpendicular_v2 = perpendicular_v2 * f1 / r_jk;
        coord_t force2 = perpendicular_v2 * dv_per_angle;
        add_force(dvelocities[ak].x, force2.x);
        add_force(dvelocities[ak].y, force2.y);
        add_force(dvelocities[ak].z, force2.z);

        add_force(dvelocities[aj].x, -(force1.x + force2.x));
        add_force(dvelocities[aj].y, -(force1.y + force2.y));
        add_force(dvelocities[aj].z, -(force1.z + force2.z));
    }
    ctx.E_bond_p.Uangle = energy_from_accum(e[bonded_p_slot()]);
    ctx.E_bond_w.Uangle = energy_from_accum(e[bonded_w_slot()]);
}

void CpuBondedForce::calc_torsions(Context& ctx) {
    /*
    Use SPFP in here.
    */
    auto* ids = data_.torsion.ids->cpu_data_p;
    auto* params = data_.torsion.params->cpu_data_p;
    auto* eslot = data_.torsion.eslot->cpu_data_p;
    auto* dvelocities = ctx.dvelocities->cpu_data_p;
    auto* coords = ctx.coords->cpu_data_p;

    energy_accum_t e[2] = {0, 0};

    for (int i = 0; i < data_.torsion.n; i++) {
        const int ai = ids[i].i, aj = ids[i].j, ak = ids[i].k, al = ids[i].l;
        real_t3 aji = real3_cast<real_t>(coords[ai] - coords[aj]);
        real_t3 ajk = real3_cast<real_t>(coords[ak] - coords[aj]);
        real_t3 akl = real3_cast<real_t>(coords[al] - coords[ak]);

        real_t3 n_ijk = cross(aji, ajk);  // The order is important
        real_t3 n_jkl = cross(akl, ajk);

        real_t norm_ijk = norm(n_ijk);
        real_t norm_jkl = norm(n_jkl);
        real_t cos_phi = dot(n_ijk, n_jkl) / (norm_ijk * norm_jkl);
        cos_phi = std::min(cos_phi, static_cast<real_t>(1.0));
        cos_phi = std::max(cos_phi, static_cast<real_t>(-1.0));
        real_t phi = acos(cos_phi);

        //
        if (dot(ajk, cross(n_ijk, n_jkl)) < 0) {
            phi = -phi;
        }

        const real_t k = params[i].k, gamma = params[i].d, paths = params[i].paths;
        const int n = static_cast<int>(params[i].n);

        auto [v, dv_per_angle] = calc_torsion(k, n, phi, gamma, paths);
        add_energy(e[eslot[i]], v);

        real_t3 di = get_perpendicular_vector(n_ijk, n_jkl);

        real_t sin_phi = sin(phi);
        real_t force_per_cos;

        if (std::abs(sin_phi) < tm06) {
            real_t cos_phi_limit = cos_phi < 0 ? -1 : 1;
            force_per_cos = k * n * n * paths * cos(n * phi - gamma) / cos_phi_limit;
        } else {
            force_per_cos = -dv_per_angle / sin_phi;
        }

        di = di / norm_ijk;

        real_t3 dl = get_perpendicular_vector(n_jkl, n_ijk);
        dl = dl / norm_jkl;

        real_t3 dpi = cross(ajk, di);
        real_t3 dpl = cross(ajk, dl);
        real_t3 dpj = cross((aji - ajk), di) + cross(akl, dl);
        real_t3 dpk = cross(static_cast<real_t>(-1) * (ajk + akl), dl) - cross(aji, di);

        add_force(dvelocities[ai].x, force_per_cos * dpi.x);
        add_force(dvelocities[ai].y, force_per_cos * dpi.y);
        add_force(dvelocities[ai].z, force_per_cos * dpi.z);

        add_force(dvelocities[al].x, force_per_cos * dpl.x);
        add_force(dvelocities[al].y, force_per_cos * dpl.y);
        add_force(dvelocities[al].z, force_per_cos * dpl.z);

        add_force(dvelocities[aj].x, force_per_cos * dpj.x);
        add_force(dvelocities[aj].y, force_per_cos * dpj.y);
        add_force(dvelocities[aj].z, force_per_cos * dpj.z);

        add_force(dvelocities[ak].x, force_per_cos * dpk.x);
        add_force(dvelocities[ak].y, force_per_cos * dpk.y);
        add_force(dvelocities[ak].z, force_per_cos * dpk.z);
    }
    ctx.E_bond_p.Utor = energy_from_accum(e[bonded_p_slot()]);
    ctx.E_bond_w.Utor = energy_from_accum(e[bonded_w_slot()]);
}

void CpuBondedForce::calc_impropers(Context& ctx) {
    /*
    The deviative part is same as torsion.

    */
    auto* ids = data_.improper.ids->cpu_data_p;
    auto* params = data_.improper.params->cpu_data_p;
    auto* eslot = data_.improper.eslot->cpu_data_p;
    auto* dvelocities = ctx.dvelocities->cpu_data_p;
    auto* coords = ctx.coords->cpu_data_p;

    energy_accum_t e[2] = {0, 0};

    for (int i = 0; i < data_.improper.n; i++) {
        const int ai = ids[i].i, aj = ids[i].j, ak = ids[i].k, al = ids[i].l;
        coord_t aji = coords[ai] - coords[aj];
        coord_t ajk = coords[ak] - coords[aj];
        coord_t akl = coords[al] - coords[ak];

        coord_t n_ijk = cross(aji, ajk);
        coord_t n_jkl = cross(akl, ajk);

        double norm_ijk = norm(n_ijk);
        double norm_jkl = norm(n_jkl);
        double cos_phi = dot(n_ijk, n_jkl) / (norm_ijk * norm_jkl);
        cos_phi = std::min(cos_phi, 1.0);
        cos_phi = std::max(cos_phi, -1.0);
        double phi = acos(cos_phi);

        if (dot(ajk, cross(n_ijk, n_jkl)) < 0) {
            phi = -phi;
        }

        const double k = params[i].x, phi0 = params[i].y;
        auto [v, dv_per_angle] = calc_improper2(k, phi0, phi);
        add_energy(e[eslot[i]], v);

        coord_t di = get_perpendicular_vector(n_ijk, n_jkl);
        double sin_phi = sin(phi);
        double force_per_cos;
        if (std::abs(sin_phi) < tm06) {
            double cos_phi_limit = cos_phi < 0 ? -1 : 1;
            force_per_cos = 4.0 * k * cos(2.0 * phi - phi0) / cos_phi_limit;
        } else {
            force_per_cos = -dv_per_angle / sin_phi;
        }
        di = di / norm_ijk;

        coord_t dl = get_perpendicular_vector(n_jkl, n_ijk);
        dl = dl / norm_jkl;

        coord_t dpi = cross(ajk, di);
        coord_t dpl = cross(ajk, dl);
        coord_t dpj = cross((aji - ajk), di) + cross(akl, dl);
        coord_t dpk = cross(-1 * (ajk + akl), dl) - cross(aji, di);

        add_force(dvelocities[ai].x, force_per_cos * dpi.x);
        add_force(dvelocities[ai].y, force_per_cos * dpi.y);
        add_force(dvelocities[ai].z, force_per_cos * dpi.z);

        add_force(dvelocities[al].x, force_per_cos * dpl.x);
        add_force(dvelocities[al].y, force_per_cos * dpl.y);
        add_force(dvelocities[al].z, force_per_cos * dpl.z);

        add_force(dvelocities[aj].x, force_per_cos * dpj.x);
        add_force(dvelocities[aj].y, force_per_cos * dpj.y);
        add_force(dvelocities[aj].z, force_per_cos * dpj.z);

        add_force(dvelocities[ak].x, force_per_cos * dpk.x);
        add_force(dvelocities[ak].y, force_per_cos * dpk.y);
        add_force(dvelocities[ak].z, force_per_cos * dpk.z);
    }

    ctx.E_bond_p.Uimp = energy_from_accum(e[bonded_p_slot()]);
    ctx.E_bond_w.Uimp = energy_from_accum(e[bonded_w_slot()]);
}