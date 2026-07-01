#include "constants.h"
#include "cuda_bonded_force.cuh"
#include "cuda_force_accumulation.cuh"
#include "geometry.h"

namespace {

__device__ void compute_bond(int i, const bond_idx_t* ids, const dparam2_t* params, const int* eslot, const coord_t* coords, dvel_t* dvel, energy_accum_t* e_bond) {
    const int ai = ids[i].i, aj = ids[i].j;
    const double k = params[i].x, r_eq = params[i].y;

    coord_t d = coords[aj] - coords[ai];
    const double r = norm(d);
    auto [v, dv] = calc_bond(k, r, r_eq);
    const double ampl = dv / r;

    atomic_add_energy(&e_bond[eslot[i]], v);
    atomic_add_force(&dvel[ai].x, -ampl * d.x);
    atomic_add_force(&dvel[ai].y, -ampl * d.y);
    atomic_add_force(&dvel[ai].z, -ampl * d.z);
    atomic_add_force(&dvel[aj].x, ampl * d.x);
    atomic_add_force(&dvel[aj].y, ampl * d.y);
    atomic_add_force(&dvel[aj].z, ampl * d.z);
}

__device__ void compute_angle(int i, const angle_idx_t* ids, const dparam2_t* params, const int* eslot, const coord_t* coords, dvel_t* dvel, energy_accum_t* e_angle) {
    const int ai = ids[i].i, aj = ids[i].j, ak = ids[i].k;
    const double k = params[i].x;
    const double th_eq = params[i].y;

    coord_t aji = coords[ai] - coords[aj];
    coord_t ajk = coords[ak] - coords[aj];

    double r_ji = norm(aji);
    double r_jk = norm(ajk);

    double cos_th = dot(aji, ajk) / (r_ji * r_jk);
    cos_th = min(cos_th, pt999);
    cos_th = max(cos_th, -pt999);

    const double th = acos(cos_th);
    auto [v, dv_per_angle] = calc_angle(k, th, th_eq);

    atomic_add_energy(&e_angle[eslot[i]], v);

    coord_t perpendicular_v = get_perpendicular_vector(aji, ajk);
    double sin_th = sin(th);
    double f1 = sin_th;
    f1 = -1.0 / f1;
    perpendicular_v = perpendicular_v * f1 / r_ji;
    coord_t force1 = perpendicular_v * dv_per_angle;
    atomic_add_force(&dvel[ai].x, force1.x);
    atomic_add_force(&dvel[ai].y, force1.y);
    atomic_add_force(&dvel[ai].z, force1.z);

    coord_t perpendicular_v2 = get_perpendicular_vector(ajk, aji);
    perpendicular_v2 = perpendicular_v2 * f1 / r_jk;
    coord_t force2 = perpendicular_v2 * dv_per_angle;
    atomic_add_force(&dvel[ak].x, force2.x);
    atomic_add_force(&dvel[ak].y, force2.y);
    atomic_add_force(&dvel[ak].z, force2.z);

    atomic_add_force(&dvel[aj].x, -(force1.x + force2.x));
    atomic_add_force(&dvel[aj].y, -(force1.y + force2.y));
    atomic_add_force(&dvel[aj].z, -(force1.z + force2.z));
}

__device__ void compute_torsion(int i, const dihe_idx_t* ids, const torsion_param_t* params, const int* eslot, const coord_t* coords, dvel_t* dvelocities, energy_accum_t* e_torsion) {
    const int ai = ids[i].i, aj = ids[i].j, ak = ids[i].k, al = ids[i].l;
    real_t3 aji = real3_cast<real_t>(coords[ai] - coords[aj]);
    real_t3 ajk = real3_cast<real_t>(coords[ak] - coords[aj]);
    real_t3 akl = real3_cast<real_t>(coords[al] - coords[ak]);

    real_t3 n_ijk = cross(aji, ajk);  // The order is important
    real_t3 n_jkl = cross(akl, ajk);

    real_t norm_ijk = norm(n_ijk);
    real_t norm_jkl = norm(n_jkl);
    real_t cos_phi = dot(n_ijk, n_jkl) / (norm_ijk * norm_jkl);
    cos_phi = min(cos_phi, static_cast<real_t>(1.0));
    cos_phi = max(cos_phi, static_cast<real_t>(-1.0));
    real_t phi = acos(cos_phi);

    //
    if (dot(ajk, cross(n_ijk, n_jkl)) < 0) {
        phi = -phi;
    }

    const real_t k = params[i].k, gamma = params[i].d, paths = params[i].paths;
    const int n = static_cast<int>(params[i].n);

    auto [v, dv_per_angle] = calc_torsion(k, n, phi, gamma, paths);
    atomic_add_energy(&e_torsion[eslot[i]], v);

    real_t3 di = get_perpendicular_vector(n_ijk, n_jkl);

    real_t sin_phi = sin(phi);

    real_t force_per_cos;
    if (abs(sin_phi) < tm06) {
        /*
        dv_per_angle = dU / dphi = -k * n * sin(n * phi - gamma) * paths
        force_per_cos = -dv_per_angle / sin_phi =  k * n * paths * sin(n * phi - gamma) / sin(phi)
        when phi = 0 or pi. sin(phi) = 0
        l'Hopital:
        f(phi) -> 0, g(phi) -> 0
        lim f(phi) / g(phi) = lim f'(phi) / g'(phi)

        lim sin(n * phi - gamma) / sin(phi)
        f(phi) = sin(n * phi - gamma)
        g(phi) = sin(phi)
        f'(phi) = n * cos(n * phi - gamma)
        g'(phi) = cos(phi)
        */

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

    atomic_add_force(&dvelocities[ai].x, force_per_cos * dpi.x);
    atomic_add_force(&dvelocities[ai].y, force_per_cos * dpi.y);
    atomic_add_force(&dvelocities[ai].z, force_per_cos * dpi.z);

    atomic_add_force(&dvelocities[al].x, force_per_cos * dpl.x);
    atomic_add_force(&dvelocities[al].y, force_per_cos * dpl.y);
    atomic_add_force(&dvelocities[al].z, force_per_cos * dpl.z);

    atomic_add_force(&dvelocities[aj].x, force_per_cos * dpj.x);
    atomic_add_force(&dvelocities[aj].y, force_per_cos * dpj.y);
    atomic_add_force(&dvelocities[aj].z, force_per_cos * dpj.z);

    atomic_add_force(&dvelocities[ak].x, force_per_cos * dpk.x);
    atomic_add_force(&dvelocities[ak].y, force_per_cos * dpk.y);
    atomic_add_force(&dvelocities[ak].z, force_per_cos * dpk.z);
}

__device__ void compute_improper(int i, const dihe_idx_t* ids, const dparam2_t* params, const int* eslot, const coord_t* coords, dvel_t* dvelocities, energy_accum_t* e_improper) {
    const int ai = ids[i].i, aj = ids[i].j, ak = ids[i].k, al = ids[i].l;
    coord_t aji = coords[ai] - coords[aj];
    coord_t ajk = coords[ak] - coords[aj];
    coord_t akl = coords[al] - coords[ak];

    coord_t n_ijk = cross(aji, ajk);
    coord_t n_jkl = cross(akl, ajk);

    double norm_ijk = norm(n_ijk);
    double norm_jkl = norm(n_jkl);
    double cos_phi = dot(n_ijk, n_jkl) / (norm_ijk * norm_jkl);
    cos_phi = min(cos_phi, 1.0);
    cos_phi = max(cos_phi, -1.0);
    double phi = acos(cos_phi);

    if (dot(ajk, cross(n_ijk, n_jkl)) < 0) {
        phi = -phi;
    }

    const double k = params[i].x, phi0 = params[i].y;
    auto [v, dv_per_angle] = calc_improper2(k, phi0, phi);
    atomic_add_energy(&e_improper[eslot[i]], v);

    coord_t di = get_perpendicular_vector(n_ijk, n_jkl);
    double sin_phi = sin(phi);
    double force_per_cos;
    if (abs(sin_phi) < tm06) {
        /*
        dv_per_angle = k * -2.0 * sin(2.0 * phi - phi0)
        force_per_cos = -dv_per_angle / sin_phi = 2.0 * k * sin(2.0 * phi - phi0) / sin(phi)

        f(phi) = 2.0 * k * sin(2.0 * phi - phi0)
        f' = 2.0 * k * cos(2.0 * phi - phi0) * 2.0

        g(phi) = sin(phi)
        g' = cos(phi)
        */
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

    atomic_add_force(&dvelocities[ai].x, force_per_cos * dpi.x);
    atomic_add_force(&dvelocities[ai].y, force_per_cos * dpi.y);
    atomic_add_force(&dvelocities[ai].z, force_per_cos * dpi.z);

    atomic_add_force(&dvelocities[al].x, force_per_cos * dpl.x);
    atomic_add_force(&dvelocities[al].y, force_per_cos * dpl.y);
    atomic_add_force(&dvelocities[al].z, force_per_cos * dpl.z);

    atomic_add_force(&dvelocities[aj].x, force_per_cos * dpj.x);
    atomic_add_force(&dvelocities[aj].y, force_per_cos * dpj.y);
    atomic_add_force(&dvelocities[aj].z, force_per_cos * dpj.z);

    atomic_add_force(&dvelocities[ak].x, force_per_cos * dpk.x);
    atomic_add_force(&dvelocities[ak].y, force_per_cos * dpk.y);
    atomic_add_force(&dvelocities[ak].z, force_per_cos * dpk.z);
}

__global__ void bonded_kernel(
    int n_bond, int n_angle, int n_torsion, int n_improper,
    const bond_idx_t* bond_ids, const dparam2_t* bond_p, const int* bond_es,
    const angle_idx_t* angl_ids, const dparam2_t* angl_p, const int* angl_es,
    const dihe_idx_t* tor_ids, const torsion_param_t* tor_p, const int* tor_es,
    const dihe_idx_t* imp_ids, const dparam2_t* imp_p, const int* imp_es,
    const coord_t* coords, dvel_t* dvel,
    energy_accum_t* e_bond, energy_accum_t* e_angle,
    energy_accum_t* e_torsion, energy_accum_t* e_improper) {
    int tid = blockIdx.x * blockDim.x + threadIdx.x;
    if (tid < n_bond) {
        compute_bond(tid, bond_ids, bond_p, bond_es, coords, dvel, e_bond);
        return;
    }
    tid -= n_bond;
    if (tid < n_angle) {
        compute_angle(tid, angl_ids, angl_p, angl_es, coords, dvel, e_angle);
        return;
    }
    tid -= n_angle;
    if (tid < n_torsion) {
        compute_torsion(tid, tor_ids, tor_p, tor_es, coords, dvel, e_torsion);
        return;
    }
    tid -= n_torsion;
    if (tid < n_improper) {
        compute_improper(tid, imp_ids, imp_p, imp_es, coords, dvel, e_improper);
        return;
    }
}
}  // namespace

void CudaBondedForce::init_backend(Context& ctx) {
    const int slots = bonded_region_slots();
    e_bond_ = std::make_unique<HostDeviceBuffer<energy_accum_t>>(slots);
    e_angle_ = std::make_unique<HostDeviceBuffer<energy_accum_t>>(slots);
    e_torsion_ = std::make_unique<HostDeviceBuffer<energy_accum_t>>(slots);
    e_improper_ = std::make_unique<HostDeviceBuffer<energy_accum_t>>(slots);
}

void CudaBondedForce::calc(Context& ctx) {
    if (!enabled()) return;

    auto zero = [](auto& buf) {
        cudaMemset(buf->gpu_data_p, 0, sizeof(energy_accum_t) * buf->length);
    };
    zero(e_bond_);
    zero(e_angle_);
    zero(e_torsion_);
    zero(e_improper_);

    const int n_bond = data_.bond.n;
    const int n_angle = data_.angle.n;
    const int n_torsion = data_.torsion.n;
    const int n_improper = data_.improper.n;
    const int total = n_bond + n_angle + n_torsion + n_improper;

    const int block = 256;
    const int grid = (total + block - 1) / block;

    bonded_kernel<<<grid, block>>>(
        n_bond, n_angle, n_torsion, n_improper,
        data_.bond.ids->gpu_data_p, data_.bond.params->gpu_data_p, data_.bond.eslot->gpu_data_p,
        data_.angle.ids->gpu_data_p, data_.angle.params->gpu_data_p, data_.angle.eslot->gpu_data_p,
        data_.torsion.ids->gpu_data_p, data_.torsion.params->gpu_data_p, data_.torsion.eslot->gpu_data_p,
        data_.improper.ids->gpu_data_p, data_.improper.params->gpu_data_p, data_.improper.eslot->gpu_data_p,
        ctx.coords->gpu_data_p, ctx.dvelocities->gpu_data_p,
        e_bond_->gpu_data_p, e_angle_->gpu_data_p, e_torsion_->gpu_data_p, e_improper_->gpu_data_p);

    e_bond_->download();
    e_angle_->download();
    e_torsion_->download();
    e_improper_->download();

    ctx.E_bond_p.Ubond = energy_from_accum(e_bond_->cpu_data_p[bonded_p_slot()]);
    ctx.E_bond_w.Ubond = energy_from_accum(e_bond_->cpu_data_p[bonded_w_slot()]);
    ctx.E_bond_p.Uangle = energy_from_accum(e_angle_->cpu_data_p[bonded_p_slot()]);
    ctx.E_bond_w.Uangle = energy_from_accum(e_angle_->cpu_data_p[bonded_w_slot()]);
    ctx.E_bond_p.Utor = energy_from_accum(e_torsion_->cpu_data_p[bonded_p_slot()]);
    ctx.E_bond_w.Utor = energy_from_accum(e_torsion_->cpu_data_p[bonded_w_slot()]);
    ctx.E_bond_p.Uimp = energy_from_accum(e_improper_->cpu_data_p[bonded_p_slot()]);
    ctx.E_bond_w.Uimp = energy_from_accum(e_improper_->cpu_data_p[bonded_w_slot()]);
}