#include <vector>

#include "cuda_settle.cuh"
#include "md_types.h"

namespace {

const int kShakeThreads = 64;
__device__ double clamp_nonnegative(double value) {
    return value > 0.0 ? value : 0.0;
}
__global__ void calc_fast_water_shake_kernel(
    int n_fast_waters,
    SettleFastWater* fast_waters,
    coord_t* coords,
    const coord_t* xcoords) {
    /*
    todo: Need to think...
    */
    int idx = blockIdx.x * blockDim.x + threadIdx.x;
    if (idx >= n_fast_waters) return;

    auto water = fast_waters[idx];
    const int o = water.o;
    const int h1 = water.h1;
    const int h2 = water.h2;

    const double xo0 = xcoords[o].x;
    const double yo0 = xcoords[o].y;
    const double zo0 = xcoords[o].z;
    const double xh10 = xcoords[h1].x;
    const double yh10 = xcoords[h1].y;
    const double zh10 = xcoords[h1].z;
    const double xh20 = xcoords[h2].x;
    const double yh20 = xcoords[h2].y;
    const double zh20 = xcoords[h2].z;

    const double xo1 = coords[o].x;
    const double yo1 = coords[o].y;
    const double zo1 = coords[o].z;
    const double xh11 = coords[h1].x;
    const double yh11 = coords[h1].y;
    const double zh11 = coords[h1].z;
    const double xh21 = coords[h2].x;
    const double yh21 = coords[h2].y;
    const double zh21 = coords[h2].z;

    const double xb0 = xh10 - xo0;
    const double yb0 = yh10 - yo0;
    const double zb0 = zh10 - zo0;
    const double xc0 = xh20 - xo0;
    const double yc0 = yh20 - yo0;
    const double zc0 = zh20 - zo0;

    const double xcom = xo1 * water.wo_div_wohh + (xh11 + xh21) * water.wh_div_wohh;
    const double ycom = yo1 * water.wo_div_wohh + (yh11 + yh21) * water.wh_div_wohh;
    const double zcom = zo1 * water.wo_div_wohh + (zh11 + zh21) * water.wh_div_wohh;

    const double xa1 = xo1 - xcom;
    const double ya1 = yo1 - ycom;
    const double za1 = zo1 - zcom;
    const double xb1 = xh11 - xcom;
    const double yb1 = yh11 - ycom;
    const double zb1 = zh11 - zcom;
    const double xc1 = xh21 - xcom;
    const double yc1 = yh21 - ycom;
    const double zc1 = zh21 - zcom;

    const double xakszd = yb0 * zc0 - zb0 * yc0;
    const double yakszd = zb0 * xc0 - xb0 * zc0;
    const double zakszd = xb0 * yc0 - yb0 * xc0;
    const double xaksxd = ya1 * zakszd - za1 * yakszd;
    const double yaksxd = za1 * xakszd - xa1 * zakszd;
    const double zaksxd = xa1 * yakszd - ya1 * xakszd;
    const double xaksyd = yakszd * zaksxd - zakszd * yaksxd;
    const double yaksyd = zakszd * xaksxd - xakszd * zaksxd;
    const double zaksyd = xakszd * yaksxd - yakszd * xaksxd;

    const double ax2 = xaksxd * xaksxd + yaksxd * yaksxd + zaksxd * zaksxd;
    const double ay2 = xaksyd * xaksyd + yaksyd * yaksyd + zaksyd * zaksyd;
    const double az2 = xakszd * xakszd + yakszd * yakszd + zakszd * zakszd;
    if (ax2 <= 0.0 || ay2 <= 0.0 || az2 <= 0.0) return;

    const double axlng_inv = rsqrt(ax2);
    const double aylng_inv = rsqrt(ay2);
    const double azlng_inv = rsqrt(az2);

    const double trns11 = xaksxd * axlng_inv;
    const double trns21 = yaksxd * axlng_inv;
    const double trns31 = zaksxd * axlng_inv;
    const double trns12 = xaksyd * aylng_inv;
    const double trns22 = yaksyd * aylng_inv;
    const double trns32 = zaksyd * aylng_inv;
    const double trns13 = xakszd * azlng_inv;
    const double trns23 = yakszd * azlng_inv;
    const double trns33 = zakszd * azlng_inv;

    const double xb0d = trns11 * xb0 + trns21 * yb0 + trns31 * zb0;
    const double yb0d = trns12 * xb0 + trns22 * yb0 + trns32 * zb0;
    const double xc0d = trns11 * xc0 + trns21 * yc0 + trns31 * zc0;
    const double yc0d = trns12 * xc0 + trns22 * yc0 + trns32 * zc0;
    const double za1d = trns13 * xa1 + trns23 * ya1 + trns33 * za1;
    const double xb1d = trns11 * xb1 + trns21 * yb1 + trns31 * zb1;
    const double yb1d = trns12 * xb1 + trns22 * yb1 + trns32 * zb1;
    const double zb1d = trns13 * xb1 + trns23 * yb1 + trns33 * zb1;
    const double xc1d = trns11 * xc1 + trns21 * yc1 + trns31 * zc1;
    const double yc1d = trns12 * xc1 + trns22 * yc1 + trns32 * zc1;
    const double zc1d = trns13 * xc1 + trns23 * yc1 + trns33 * zc1;

    const double sinphi = za1d * water.ra_inv;
    const double cosphi = sqrt(clamp_nonnegative(1.0 - sinphi * sinphi));
    if (cosphi <= 0.0) return;
    const double sinpsi = (zb1d - zc1d) / (water.rhh * cosphi);
    const double cospsi = sqrt(clamp_nonnegative(1.0 - sinpsi * sinpsi));

    const double ya2d = water.ra * cosphi;
    double xb2d = -water.rc * cospsi;
    const double yb2d = -water.rb * cosphi - water.rc * sinpsi * sinphi;
    const double yc2d = -water.rb * cosphi + water.rc * sinpsi * sinphi;
    xb2d = -0.5 *
           sqrt(clamp_nonnegative(water.rhh2 - (yb2d - yc2d) * (yb2d - yc2d) - (zb1d - zc1d) * (zb1d - zc1d)));

    const double alpa = xb2d * (xb0d - xc0d) + yb0d * yb2d + yc0d * yc2d;
    const double beta = xb2d * (yc0d - yb0d) + xb0d * yb2d + xc0d * yc2d;
    const double gama = xb0d * yb1d - xb1d * yb0d + xc0d * yc1d - xc1d * yc0d;
    const double al2be2 = alpa * alpa + beta * beta;
    if (al2be2 <= 0.0) return;

    const double sinthe = (alpa * gama - beta * sqrt(clamp_nonnegative(al2be2 - gama * gama))) / al2be2;
    const double costhe = sqrt(clamp_nonnegative(1.0 - sinthe * sinthe));

    const double xa3d = -ya2d * sinthe;
    const double ya3d = ya2d * costhe;
    const double za3d = za1d;
    const double xb3d = xb2d * costhe - yb2d * sinthe;
    const double yb3d = xb2d * sinthe + yb2d * costhe;
    const double zb3d = zb1d;
    const double xc3d = -xb2d * costhe - yc2d * sinthe;
    const double yc3d = -xb2d * sinthe + yc2d * costhe;
    const double zc3d = zc1d;

    coords[o].x = xcom + trns11 * xa3d + trns12 * ya3d + trns13 * za3d;
    coords[o].y = ycom + trns21 * xa3d + trns22 * ya3d + trns23 * za3d;
    coords[o].z = zcom + trns31 * xa3d + trns32 * ya3d + trns33 * za3d;
    coords[h1].x = xcom + trns11 * xb3d + trns12 * yb3d + trns13 * zb3d;
    coords[h1].y = ycom + trns21 * xb3d + trns22 * yb3d + trns23 * zb3d;
    coords[h1].z = zcom + trns31 * xb3d + trns32 * yb3d + trns33 * zb3d;
    coords[h2].x = xcom + trns11 * xc3d + trns12 * yc3d + trns13 * zc3d;
    coords[h2].y = ycom + trns21 * xc3d + trns22 * yc3d + trns23 * zc3d;
    coords[h2].z = zcom + trns31 * xc3d + trns32 * yc3d + trns33 * zc3d;
}

}  // namespace

SettleFastWater make_fast_water(int o, int h1, int h2, double roh2, double rhh2, double wo, double wh) {
    const double roh = sqrt(roh2);
    const double rhh = sqrt(rhh2);
    const double height = sqrt(roh2 - 0.25 * rhh2);
    const double total_mass = wo + 2.0 * wh;
    const double com_y = 2.0 * wh * height / total_mass;

    return SettleFastWater{
        o, h1, h2, com_y, 1.0 / com_y, height - com_y, 0.5 * rhh, rhh, rhh * rhh, wo / total_mass, wh / total_mass};
}

void CudaSettleSolver::init(const std::vector<SettleFastWater>& settle_fast_waters) {
    settle_fast_waters_ = HostDeviceBuffer<SettleFastWater>::from_vector(settle_fast_waters, true);
}

void CudaSettleSolver::apply(coord_t* d_coords, const coord_t* d_xcoords, const double* d_winv) {
    if (settle_fast_waters_->length > 0) {
        const int grid_blocks = (settle_fast_waters_->length + kShakeThreads - 1) / kShakeThreads;
        calc_fast_water_shake_kernel<<<grid_blocks, kShakeThreads>>>(settle_fast_waters_->length, settle_fast_waters_->gpu_data_p, d_coords, d_xcoords);
    }
}

void CudaSettle::init_backend(Context& ctx) {
    /*
    Settle can only be used in water.
    */
    int bond_num = data_.n_constraints;
    std::vector<ConstraintBond> constraint_bonds(data_.constraint_bonds->cpu_data_p, data_.constraint_bonds->cpu_data_p + bond_num);
    init_from_bonds(ctx, constraint_bonds);
}

void CudaSettle::init_from_bonds(Context& ctx, const std::vector<ConstraintBond>& bonds) {
    std::vector<SettleFastWater> settle_fast_waters;
    const auto* winv = ctx.winv->cpu_data_p;

    std::map<std::pair<int, int>, int> bonds_lookup_table;
    for (int i = 0; i < bonds.size(); i++) {
        int ai = bonds[i].ai - 1;
        int aj = bonds[i].aj - 1;
        bonds_lookup_table[std::minmax(ai, aj)] = i;
    }

    auto find_bond_index = [&](int ai, int aj) {
        auto it = bonds_lookup_table.find(std::minmax(ai, aj));
        return it == bonds_lookup_table.end() ? -1 : it->second;
    };

    for (int w = 0; w < ctx.n_waters(); w++) {
        const int o = ctx.n_atoms_solute + 3 * w;
        const int h1 = o + 1;
        const int h2 = o + 2;

        const int oh1 = find_bond_index(o, h1);
        const int oh2 = find_bond_index(o, h2);
        const int hh = find_bond_index(h1, h2);

        if (oh1 < 0 || oh2 < 0 || hh < 0) continue;
        if (bonds[oh1].dist2 != bonds[oh2].dist2) {
            printf("[Settle] two o-h bonds dists are not equal. Wrong config!");
            std::exit(EXIT_FAILURE);
        }
        if (bonds[oh1].dist2 <= 0.25 * bonds[hh].dist2) {
            // triangle height from oxygen to the H-H midpoint height = sqrt(rOH * rOH - 0.25 * rHH * rHH);
            printf("[Settle] Triangle height from oxygen to the H-H midpoint height is smaller than 0. Wrong config");
            std::exit(EXIT_FAILURE);
        }

        const double wo = 1.0 / winv[o];
        const double wh1 = 1.0 / winv[h1];
        const double wh2 = 1.0 / winv[h2];
        if (wh1 != wh2) {
            continue;
        }
        settle_fast_waters.push_back(make_fast_water(o, h1, h2, bonds[oh1].dist2, bonds[hh].dist2, wo, wh1));
    }
    solver_.init(settle_fast_waters);
}

void CudaSettle::apply(Context& ctx, HostDeviceBuffer<coord_t>& xcoords) {
    apply_to(ctx, ctx.coords->gpu_data_p, xcoords.gpu_data_p);
}

void CudaSettle::apply_to(Context& ctx, coord_t* d_coords, coord_t* d_xcoords) {
    solver_.apply(d_coords, d_xcoords, ctx.winv->gpu_data_p);
}

void CudaSettle::initial_constraint(Context& ctx) {
    HostDeviceBuffer<coord_t> xcoords(ctx.n_atoms, true, true);
    coord_t* d_coords = ctx.coords->gpu_data_p;
    coord_t* d_xcoords = xcoords.gpu_data_p;

    check_cuda(cudaMemcpy(
        d_xcoords,
        d_coords,
        ctx.n_atoms * sizeof(coord_t),
        cudaMemcpyDeviceToDevice));
    apply_to(ctx, d_coords, d_xcoords);

    ctx.coords->download();
    coord_t* coords = ctx.coords->cpu_data_p;
    vel_t* velocities = ctx.velocities->cpu_data_p;
    coord_t* host_xcoords = xcoords.cpu_data_p;
    for (int i = 0; i < ctx.n_atoms; i++) {
        host_xcoords[i].x = coords[i].x - ctx.dt * velocities[i].x;
        host_xcoords[i].y = coords[i].y - ctx.dt * velocities[i].y;
        host_xcoords[i].z = coords[i].z - ctx.dt * velocities[i].z;
    }

    xcoords.upload();
    apply_to(ctx, d_xcoords, d_coords);
    xcoords.download();

    for (int i = 0; i < ctx.n_atoms; i++) {
        velocities[i].x = (coords[i].x - host_xcoords[i].x) / ctx.dt;
        velocities[i].y = (coords[i].y - host_xcoords[i].y) / ctx.dt;
        velocities[i].z = (coords[i].z - host_xcoords[i].z) / ctx.dt;
    }
    ctx.velocities->upload();
}
