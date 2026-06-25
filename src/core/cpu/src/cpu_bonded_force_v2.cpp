#include "cpu_bonded_force_v2.h"

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
        cos_th = std::min(cos_th, 1.0);
        cos_th = std::max(cos_th, -1.0);
        const double th = acos(cos_th);
        auto [v, dv_per_angle] = calc_angle(k, th, th_eq);

        add_energy(e[eslot[i]], v);

        // Now we should get a vector that is perpendicular to aji inside the plane formed by atoms i-j-k
        coord_t perpendicular_v = get_perpendicular_vector(aji, ajk);
        perpendicular_v = perpendicular_v / std::max(norm(perpendicular_v), k_singular_sin_epsilon);  // unit, guarded

        // get_perpendicular_vector points toward k, but increasing theta moves i away from k, so flip it
        perpendicular_v = perpendicular_v * -1;
        // Change energy per angle to energy per distance.
        // Distance = R * theta (arc length = radius * angle)
        double dv = dv_per_angle / r_ji;
        coord_t force1 = perpendicular_v * dv;
        add_force(dvelocities[ai].x, force1.x);
        add_force(dvelocities[ai].y, force1.y);
        add_force(dvelocities[ai].z, force1.z);

        coord_t perpendicular_v2 = get_perpendicular_vector(ajk, aji);
        perpendicular_v2 = perpendicular_v2 / std::max(norm(perpendicular_v2), k_singular_sin_epsilon);  // unit, guarded
        perpendicular_v2 = perpendicular_v2 * -1;
        dv = dv_per_angle / r_jk;
        coord_t force2 = perpendicular_v2 * dv;
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
    auto* ids = data_.torsion.ids->cpu_data_p;
    auto* params = data_.torsion.params->cpu_data_p;
    auto* eslot = data_.torsion.eslot->cpu_data_p;
    auto* dvelocities = ctx.dvelocities->cpu_data_p;
    auto* coords = ctx.coords->cpu_data_p;

    energy_accum_t e[2] = {0, 0};

    for (int i = 0; i < data_.torsion.n; i++) {
        const int ai = ids[i].i, aj = ids[i].j, ak = ids[i].k, al = ids[i].l;
        real_t3 aji = coords[ai] - coords[aj];
        real_t3 ajk = coords[ak] - coords[aj];
        real_t3 akl = coords[al] - coords[ak];

        real_t3 n_ijk = cross(aji, ajk);  // The order is important
        real_t3 n_jkl = cross(akl, ajk);

        n_ijk = n_ijk / norm(n_ijk);
        n_jkl = n_jkl / norm(n_jkl);
        real_t cos_phi = dot(n_ijk, n_jkl);
        cos_phi = std::min(cos_phi, static_cast<real_t>(1.0));
        cos_phi = std::max(cos_phi, static_cast<real_t>(-1.0));
        real_t phi = acos(cos_phi);

        const real_t k = params[i].k, gamma = params[i].d, paths = params[i].paths;
        const int n = params[i].n;

        auto [v, dv_per_angle] = calc_torsion(k, n, phi, gamma, paths);
    }
}