#include "bonded_force.h"  // calc_angle (HD)
#include "constants.h"
#include "context.h"
#include "cuda_force_accumulation.cuh"
#include "cuda_restraint_force.cuh"
#include "geometry.h"
#include "q_math.h"

namespace {

// pshell + fixed atoms
__device__ void compute_posrestr(int i, const int* ids, const double* kk,
                                 const int* eslot_a, const int* eslot_b,
                                 const coord_t* coords, const coord_t* coords_init,
                                 dvel_t* dvel, energy_accum_t* e) {
    const int idx = ids[i];
    const double k = kk[i];
    const double dx = coords[idx].x - coords_init[idx].x;
    const double dy = coords[idx].y - coords_init[idx].y;
    const double dz = coords[idx].z - coords_init[idx].z;
    const double ener = 0.5 * k * (dx * dx + dy * dy + dz * dz);

    if (eslot_a[i] >= 0) atomic_add_energy(&e[eslot_a[i]], ener);
    if (eslot_b[i] >= 0) atomic_add_energy(&e[eslot_b[i]], ener);

    atomic_add_force(&dvel[idx].x, k * dx);
    atomic_add_force(&dvel[idx].y, k * dy);
    atomic_add_force(&dvel[idx].z, k * dz);
}

// sequence restraints: one thread owns the whole [ai,aj) group
__device__ void compute_restrseq(int s, const restrseq_t* recs,
                                 const atype_t* atypes, const catype_t* catypes,
                                 const bool* heavy, const coord_t* coords,
                                 const coord_t* coords_init, dvel_t* dvel, energy_accum_t* e) {
    const double k = recs[s].k;
    coord_t d = {0, 0, 0};

    if (recs[s].to_center == 1) {
        int n_ctr = 0;
        for (int i = recs[s].ai - 1; i < recs[s].aj - 1; i++)
            if (heavy[i] || recs[s].ih) {
                n_ctr++;
                d = d + coords[i] - coords_init[i];
            }
        if (n_ctr > 0) {
            d = d * (1.0 / static_cast<double>(n_ctr));
            atomic_add_energy(&e[E_RESTR_PRES], 0.5 * k * norm2(d));
            for (int i = recs[s].ai - 1; i < recs[s].aj - 1; i++)
                if (heavy[i] || recs[s].ih) {
                    const double tmp = catypes[atypes[i].code - 1].m / CARBON_MASS;
                    atomic_add_force(&dvel[i].x, k * d.x * tmp);
                    atomic_add_force(&dvel[i].y, k * d.y * tmp);
                    atomic_add_force(&dvel[i].z, k * d.z * tmp);
                }
        }
    } else if (recs[s].to_center == 2) {
        double totmass = 0;
        for (int i = recs[s].ai - 1; i < recs[s].aj - 1; i++)
            if (heavy[i] || recs[s].ih) {
                const double mass = catypes[atypes[i].code - 1].m;
                totmass += mass;
                d = d + (coords[i] - coords_init[i]) * mass;
            }
        if (totmass > 0) {
            d = d * (1.0 / totmass);
            atomic_add_energy(&e[E_RESTR_PRES], 0.5 * k * norm2(d));
            for (int i = recs[s].ai - 1; i < recs[s].aj - 1; i++)
                if (heavy[i] || recs[s].ih) {
                    atomic_add_force(&dvel[i].x, k * d.x);
                    atomic_add_force(&dvel[i].y, k * d.y);
                    atomic_add_force(&dvel[i].z, k * d.z);
                }
        }
    } else {
        for (int i = recs[s].ai - 1; i < recs[s].aj - 1; i++)
            if (heavy[i] || recs[s].ih) {
                d = coords[i] - coords_init[i];
                atomic_add_energy(&e[E_RESTR_PRES], 0.5 * k * norm2(d));
                atomic_add_force(&dvel[i].x, k * d.x);
                atomic_add_force(&dvel[i].y, k * d.y);
                atomic_add_force(&dvel[i].z, k * d.z);
            }
    }
}

// atom restraints (FEP)
__device__ void compute_restrpos(int ir, const restrpos_t* recs, const coord_t* coords,
                                 const double* lambdas, int n_lambdas,
                                 dvel_t* dvel, energy_accum_t* e) {
    const int state = recs[ir].ipsi - 1;
    const int i = recs[ir].a - 1;
    coord_t d = coords[i] - recs[ir].x;
    const double lambda = (recs[ir].ipsi != 0) ? lambdas[state] : 1.0;
    coord_t k = recs[ir].k;
    const double ener = 0.5 * (k.x * d.x * d.x + k.y * d.y * d.y + k.z * d.z * d.z);

    atomic_add_force(&dvel[i].x, k.x * d.x * lambda);
    atomic_add_force(&dvel[i].y, k.y * d.y * lambda);
    atomic_add_force(&dvel[i].z, k.z * d.z * lambda);

    if (recs[ir].ipsi == 0) {
        for (int st = 0; st < n_lambdas; st++)
            atomic_add_energy(&e[EnergyBuffer::eq_index(ENERGY_FIXED_COUNT, st, EQ_RESTR_URESTR)], ener);
        if (n_lambdas == 0) atomic_add_energy(&e[E_RESTR_PRES], ener);
    } else {
        atomic_add_energy(&e[EnergyBuffer::eq_index(ENERGY_FIXED_COUNT, state, EQ_RESTR_URESTR)], ener);
    }
}

// distance restraints (FEP), flat-bottom [d1,d2]
__device__ void compute_restrdis(int ir, const restrdis_t* recs, const coord_t* coords,
                                 const double* lambdas, int n_lambdas,
                                 dvel_t* dvel, energy_accum_t* e) {
    const int state = recs[ir].ipsi - 1;
    const int i = recs[ir].ai - 1;
    const int j = recs[ir].aj - 1;
    coord_t aij = coords[j] - coords[i];
    const double lambda = (recs[ir].ipsi != 0) ? lambdas[state] : 1.0;
    const double r = norm(aij);

    double dr;
    if (r < recs[ir].d1)
        dr = r - recs[ir].d1;
    else if (r > recs[ir].d2)
        dr = r - recs[ir].d2;
    else
        return;  // inside flat bottom: no force/energy

    const double ener = 0.5 * recs[ir].k * dr * dr;
    const double dv = lambda * recs[ir].k * dr / r;

    atomic_add_force(&dvel[j].x, aij.x * dv);
    atomic_add_force(&dvel[j].y, aij.y * dv);
    atomic_add_force(&dvel[j].z, aij.z * dv);
    atomic_add_force(&dvel[i].x, -aij.x * dv);
    atomic_add_force(&dvel[i].y, -aij.y * dv);
    atomic_add_force(&dvel[i].z, -aij.z * dv);

    if (recs[ir].ipsi == 0) {
        for (int st = 0; st < n_lambdas; st++)
            atomic_add_energy(&e[EnergyBuffer::eq_index(ENERGY_FIXED_COUNT, st, EQ_RESTR_URESTR)], ener);
        if (n_lambdas == 0) atomic_add_energy(&e[E_RESTR_PRES], ener);
    } else {
        atomic_add_energy(&e[EnergyBuffer::eq_index(ENERGY_FIXED_COUNT, state, EQ_RESTR_URESTR)], ener);
    }
}

// angle restraints (FEP)
__device__ void compute_restrang(int ir, const restrang_t* recs, const coord_t* coords,
                                 const double* lambdas, int n_lambdas,
                                 dvel_t* dvel, energy_accum_t* e) {
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
    cos_th = min(cos_th, pt999);
    cos_th = max(cos_th, -pt999);
    const double th = acos(cos_th);

    auto res = calc_angle(recs[ir].k, th, deg2rad(recs[ir].ang));
    const double v = res.x, dv_per_angle = res.y;

    double f1 = -1.0 / sin(th);
    coord_t p1 = get_perpendicular_vector(aji, ajk) * f1 / r_ji;
    coord_t force1 = p1 * dv_per_angle * lambda;
    atomic_add_force(&dvel[i].x, force1.x);
    atomic_add_force(&dvel[i].y, force1.y);
    atomic_add_force(&dvel[i].z, force1.z);

    coord_t p2 = get_perpendicular_vector(ajk, aji) * f1 / r_jk;
    coord_t force2 = p2 * dv_per_angle * lambda;
    atomic_add_force(&dvel[k].x, force2.x);
    atomic_add_force(&dvel[k].y, force2.y);
    atomic_add_force(&dvel[k].z, force2.z);

    atomic_add_force(&dvel[j].x, -(force1.x + force2.x));
    atomic_add_force(&dvel[j].y, -(force1.y + force2.y));
    atomic_add_force(&dvel[j].z, -(force1.z + force2.z));

    if (recs[ir].ipsi == 0) {
        for (int st = 0; st < n_lambdas; st++)
            atomic_add_energy(&e[EnergyBuffer::eq_index(ENERGY_FIXED_COUNT, st, EQ_RESTR_URESTR)], v);
        if (n_lambdas == 0) atomic_add_energy(&e[E_RESTR_PRES], v);
    } else {
        atomic_add_energy(&e[EnergyBuffer::eq_index(ENERGY_FIXED_COUNT, state, EQ_RESTR_URESTR)], v);
    }
}

// wall restraints: one thread owns the whole [ai,aj) group
__device__ void compute_restrwall(int ir, const restrwall_t* recs, const bool* heavy,
                                  const coord_t* coords, coord_t solvent_center,
                                  dvel_t* dvel, energy_accum_t* e) {
    const double k = recs[ir].k;
    for (int i = recs[ir].ai - 1; i < recs[ir].aj - 1; i++) {
        if (!(heavy[i] || recs[ir].ih)) continue;
        coord_t d = coords[i] - solvent_center;
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
        atomic_add_energy(&e[EQ_RESTR_URESTR], ener);  // NOTE: preserved from CPU (flat index)
        atomic_add_force(&dvel[i].x, dv * d.x);
        atomic_add_force(&dvel[i].y, dv * d.y);
        atomic_add_force(&dvel[i].z, dv * d.z);
    }
}

__global__ void restraint_kernel(
    int n_posrestr, int n_restrseq, int n_restrpos, int n_restrdis, int n_restrang, int n_restrwall, int n_atoms, int energy_stride,
    const int* pr_ids, const double* pr_k, const int* pr_esa, const int* pr_esb,
    const restrseq_t* seq_recs, const restrpos_t* pos_recs, const restrdis_t* dis_recs,
    const restrang_t* ang_recs, const restrwall_t* wall_recs,
    const atype_t* atypes, const catype_t* catypes, const bool* heavy,
    const double* lambdas, int n_lambdas, coord_t solvent_center,
    const coord_t* coords, const coord_t* coords_init,
    dvel_t* dvel, energy_accum_t* e) {
    const int replica = blockIdx.y;
    const int atom_offset = replica * n_atoms;
    coords += atom_offset;
    dvel += atom_offset;
    e += replica * energy_stride;

    int tid = blockIdx.x * blockDim.x + threadIdx.x;
    if (tid < n_posrestr) {
        compute_posrestr(tid, pr_ids, pr_k, pr_esa, pr_esb, coords, coords_init, dvel, e);
        return;
    }
    tid -= n_posrestr;
    if (tid < n_restrseq) {
        compute_restrseq(tid, seq_recs, atypes, catypes, heavy, coords, coords_init, dvel, e);
        return;
    }
    tid -= n_restrseq;
    if (tid < n_restrpos) {
        compute_restrpos(tid, pos_recs, coords, lambdas, n_lambdas, dvel, e);
        return;
    }
    tid -= n_restrpos;
    if (tid < n_restrdis) {
        compute_restrdis(tid, dis_recs, coords, lambdas, n_lambdas, dvel, e);
        return;
    }
    tid -= n_restrdis;
    if (tid < n_restrang) {
        compute_restrang(tid, ang_recs, coords, lambdas, n_lambdas, dvel, e);
        return;
    }
    tid -= n_restrang;
    if (tid < n_restrwall) {
        compute_restrwall(tid, wall_recs, heavy, coords, solvent_center, dvel, e);
        return;
    }
}

}  // namespace

void CudaRestraintForce::init_backend(Context& ctx) {}

void CudaRestraintForce::calc(Context& ctx) {
    if (!enabled()) return;

    const int n_posrestr = data_.posrestr.n;
    const int n_restrseq = data_.restrseq.n;
    const int n_restrpos = data_.restrpos.n;
    const int n_restrdis = data_.restrdis.n;
    const int n_restrang = data_.restrang.n;
    const int n_restrwall = data_.restrwall.n;
    const int total = n_posrestr + n_restrseq + n_restrpos + n_restrdis + n_restrang + n_restrwall;

    const int block = 256;
    const int grid = (total + block - 1) / block;

    restraint_kernel<<<dim3(grid, ctx.n_replicates()), block>>>(
        n_posrestr, n_restrseq, n_restrpos, n_restrdis, n_restrang, n_restrwall, ctx.n_atoms, ctx.energy.replica_stride(),
        data_.posrestr.ids->gpu_data_p, data_.posrestr.k->gpu_data_p,
        data_.posrestr.eslot_a->gpu_data_p, data_.posrestr.eslot_b->gpu_data_p,
        n_restrseq ? data_.restrseq.recs->gpu_data_p : nullptr,
        n_restrpos ? data_.restrpos.recs->gpu_data_p : nullptr,
        n_restrdis ? data_.restrdis.recs->gpu_data_p : nullptr,
        n_restrang ? data_.restrang.recs->gpu_data_p : nullptr,
        n_restrwall ? data_.restrwall.recs->gpu_data_p : nullptr,
        ctx.atypes->gpu_data_p, ctx.catypes->gpu_data_p, ctx.heavy->gpu_data_p,
        ctx.lambdas->gpu_data_p, ctx.n_lambdas(), ctx.topo.solvent_center,
        ctx.coords->gpu_data_p, ctx.coords_init->gpu_data_p,
        ctx.dvelocities->gpu_data_p, ctx.energy.device());
}