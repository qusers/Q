#include "cuda/include/CudaContext.cuh"
#include "cuda/include/CudaNonbondedPWForce.cuh"

namespace CudaNonbondedPWForce {
// Declare any necessary static variables or device pointers here
bool is_initialized = false;
struct calc_pw_t {
    dvel_t P;
    dvel_t W;
};

double *D_PW_evdw_TOT, *D_PW_ecoul_TOT;


__device__ __forceinline__ void idx2xy(int n, int t, int& x, int& y) {
    x = (int)floorf((2 * n + 1 - sqrtf((2 * n + 1) * (2 * n + 1) - 8 * t)) * 0.5f);
    y = t - (x * n - (x * (x - 1) >> 1));
    if (y < 0) {
        x--;
        y += (n - x);
    }
    y += x;
}

__device__ __forceinline__ double shfl(double v, int srcLane, unsigned mask = 0xffffffffu) {
    int2 a = *reinterpret_cast<int2*>(&v);
    a.x = __shfl_sync(mask, a.x, srcLane);
    a.y = __shfl_sync(mask, a.y, srcLane);
    return *reinterpret_cast<double*>(&a);
}

__device__ __forceinline__ coord_t shfl_coord(coord_t v, int srcLane, unsigned mask = 0xffffffffu) {
    v.x = shfl(v.x, srcLane, mask);
    v.y = shfl(v.y, srcLane, mask);
    v.z = shfl(v.z, srcLane, mask);
    return v;
}

__device__ __forceinline__ catype_t shfl_catype(catype_t v, int srcLane, unsigned mask = 0xffffffffu) {
    v.code = __shfl_sync(mask, v.code, srcLane);
    v.m = shfl(v.m, srcLane, mask);
    v.aii_normal = shfl(v.aii_normal, srcLane, mask);
    v.bii_normal = shfl(v.bii_normal, srcLane, mask);
    v.aii_polar = shfl(v.aii_polar, srcLane, mask);
    v.bii_polar = shfl(v.bii_polar, srcLane, mask);
    v.aii_1_4 = shfl(v.aii_1_4, srcLane, mask);
    v.bii_1_4 = shfl(v.bii_1_4, srcLane, mask);
    return v;
}

__device__ void calculate_unforce_bound(
    // const int x_idx,
    // const int y_idx,

    const coord_t& x,
    const coord_t& y,

    const double x_charge,
    const double y_charge,

    const double x_aii,
    const double y_aii,

    const double x_bii,
    const double y_bii,

    const double coulomb_constant,

    const double scaling,

    double& evdw,
    double& ecoul,
    double& dv) {
    double3 d = {x.x - y.x, x.y - y.y, x.z - y.z};

    double r2 = 1.0 / (d.x * d.x + d.y * d.y + d.z * d.z);
    double r = sqrt(r2);
    double r6 = r2 * r2 * r2;

    ecoul = scaling * coulomb_constant * x_charge * y_charge * r;

    double v_a = r6 * r6 * x_aii * y_aii;
    double v_b = r6 * x_bii * y_bii;
    evdw = v_a - v_b;
    dv = r2 * (-ecoul - 12.0 * v_a + 6.0 * v_b);
}

__global__ void calc_pw(
    const int NX,
    const int NY,

    const charge_t* D_charges,
    const ccharge_t* D_ccharges,  // Charge information

    const atype_t* D_atypes,
    const catype_t* D_catypes,  // Atom type information

    const topo_t D_topo,  // Only needs coulomb_constant

    const bool* D_excluded,  // Exclusion information

    const int* D_LJ_matrix,  // bonding information, bond23 or bond14, for filtering and scaling

    const p_atom_t* D_patoms,  // the idx for X
    // Here needs the idx for Y. But W does not have such idx.

    const coord_t* coords,  // Coordinates infomation

    dvel_t* dvelocities,  // Output dvelocities

    // Other helper info, now we need it, but maybe we can optimize it later
    const int n_atoms_solute,
    bool /*x_is_solute*/,

    double* Evdw_TOT,
    double* Ecoul_TOT

) {
    const int x_blocks_num = (NX + 31) >> 5;  // solute blocks
    const int y_blocks_num = (NY + 31) >> 5;  // water blocks
    const int total_tiles = x_blocks_num * y_blocks_num;  // full rectangle
    
    const int warpsPerBlock = blockDim.x >> 5;
    const int tid = threadIdx.x;
    const int lane = tid & 31;
    const int warpInBlock = tid >> 5;

    const int tile = blockIdx.x * warpsPerBlock + warpInBlock;

    if (tile >= total_tiles) {
        return;  // extra warps in last block
    }

    int tile_y = tile / x_blocks_num;
    int tile_x = tile - tile_y * x_blocks_num;

    const int base_x = tile_x << 5;
    const int base_y = tile_y << 5;

    int x_idx = base_x + lane;
    int y_idx = base_y + lane;

    int x_atom = (x_idx < NX) ? (D_patoms[x_idx].a - 1) : -1;  // solute atom index
    int y_atom = (y_idx < NY) ? (n_atoms_solute + y_idx) : -1;  // water atom index

    coord_t invalid = {-1e9, -1e9, -1e9};
    coord_t x_coord = (x_atom >= 0) ? coords[x_atom] : invalid;
    coord_t y_coord = (y_atom >= 0) ? coords[y_atom] : invalid;

    bool x_excluded = (x_atom >= 0) ? D_excluded[x_atom] : true;
    bool y_excluded = (y_atom >= 0) ? D_excluded[y_atom] : true;

    double x_charge = (x_atom >= 0) ? D_ccharges[D_charges[x_atom].code - 1].charge : 0.0;
    double y_charge = (y_atom >= 0) ? D_ccharges[D_charges[y_atom].code - 1].charge : 0.0;  //  TODO: FIX THIS!!! WILL NOT WORK WITH QATOMS!!!!!,
                                                                                            //  For solute atom, That is fine. But for solvent atom, it is wrong.

    catype_t x_type = (x_atom >= 0) ? D_catypes[D_atypes[x_atom].code - 1] : catype_t{};
    catype_t y_type = (y_atom >= 0) ? D_catypes[D_atypes[y_atom].code - 1] : catype_t{};  // todo: Is that right for solvent?

    double3 x_force = {0.0, 0.0, 0.0};
    double3 y_force = {0.0, 0.0, 0.0};

    double evdw_sum = 0;
    double ecoul_sum = 0;

    const unsigned mask = 0xffffffffu;

    auto is_valid = [&]() -> bool {
        if (x_idx >= NX || y_idx >= NY) return false;
        if (x_excluded || y_excluded) return false;
        bool bond23 = false;
        if (bond23) return false;
        return true;
    };

    auto do_shuffle = [&]() {
        const int src = (lane + 1) & 31;
        y_idx = __shfl_sync(mask, y_idx, src);
        y_atom = __shfl_sync(mask, y_atom, src);
        y_coord = shfl_coord(y_coord, src, mask);
        y_excluded = __shfl_sync(mask, y_excluded, src);
        y_charge = __shfl_sync(mask, y_charge, src);
        y_type = shfl_catype(y_type, src, mask);

        y_force.x = __shfl_sync(mask, y_force.x, src);
        y_force.y = __shfl_sync(mask, y_force.y, src);
        y_force.z = __shfl_sync(mask, y_force.z, src);
    };

    for (int i = 0; i < 32; i++) {
        if (is_valid()) {
            // Here don't need to check LJ, LJ is only for solute-solute interactions
            bool bond14 = false;

            double scaling = bond14 ? D_topo.el14_scale : 1.0;
            double ai_aii = bond14 ? x_type.aii_1_4 : x_type.aii_normal;
            double aj_aii = bond14 ? y_type.aii_1_4 : y_type.aii_normal;

            double ai_bii = bond14 ? x_type.bii_1_4 : x_type.bii_normal;
            double aj_bii = bond14 ? y_type.bii_1_4 : y_type.bii_normal;

            double evdw = 0, ecoul = 0, dv = 0;

            calculate_unforce_bound(
                x_coord,
                y_coord,
                x_charge,
                y_charge,
                ai_aii,
                aj_aii,
                ai_bii,
                aj_bii,
                D_topo.coulomb_constant,
                scaling,
                evdw,
                ecoul,
                dv);

            evdw_sum += evdw;
            ecoul_sum += ecoul;

            double3 d = {x_coord.x - y_coord.x, x_coord.y - y_coord.y, x_coord.z - y_coord.z};
            y_force.x -= dv * d.x;
            y_force.y -= dv * d.y;
            y_force.z -= dv * d.z;

            x_force.x += dv * d.x;
            x_force.y += dv * d.y;
            x_force.z += dv * d.z;
        }
        do_shuffle();
    }

    if (x_atom >= 0) {
        atomicAdd(&dvelocities[x_atom].x, x_force.x);
        atomicAdd(&dvelocities[x_atom].y, x_force.y);
        atomicAdd(&dvelocities[x_atom].z, x_force.z);
    }
    if (y_atom >= 0) {
        atomicAdd(&dvelocities[y_atom].x, y_force.x);
        atomicAdd(&dvelocities[y_atom].y, y_force.y);
        atomicAdd(&dvelocities[y_atom].z, y_force.z);
    }

    for (int offset = 16; offset > 0; offset >>= 1) {
        evdw_sum += __shfl_down_sync(0xffffffffu, evdw_sum, offset);
        ecoul_sum += __shfl_down_sync(0xffffffffu, ecoul_sum, offset);
    }
    if (lane == 0) {
        // printf("Tile x=%d y=%d evdw_sum=%f ecoul_sum=%f\n", tile_x, tile_y, evdw_sum, ecoul_sum);
        atomicAdd(Evdw_TOT, evdw_sum);
        atomicAdd(Ecoul_TOT, ecoul_sum);
    }
}

}  // namespace CudaNonbondedPWForce

void calc_nonbonded_pw_forces_host_v2() {
    using namespace CudaNonbondedPWForce;
    int mem_size_W = 3 * n_waters * sizeof(coord_t);
    int mem_size_DV_W = 3 * n_waters * sizeof(dvel_t);
    int mem_size_DV_X = n_atoms_solute * sizeof(dvel_t);
    int mem_size_PW_MAT = 3 * n_waters * n_patoms * sizeof(calc_pw_t);

    int n_blocks_p = (n_patoms + BLOCK_SIZE - 1) / BLOCK_SIZE;
    int n_blocks_w = (3 * n_waters + BLOCK_SIZE - 1) / BLOCK_SIZE;

    int mem_size_PW_Evdw = n_blocks_p * n_blocks_w * sizeof(double);
    int mem_size_PW_Ecoul = n_blocks_p * n_blocks_w * sizeof(double);

    CudaContext& ctx = CudaContext::instance();
    auto X = ctx.d_coords;
    auto DV_X = ctx.d_dvelocities;
    auto D_ccharges = ctx.d_ccharges;
    auto D_charges = ctx.d_charges;
    auto D_catypes = ctx.d_catypes;
    auto D_atypes = ctx.d_atypes;
    auto D_patoms = ctx.d_p_atoms;
    auto D_excluded = ctx.d_excluded;

    cudaMemset(D_PW_evdw_TOT, 0, sizeof(double));
    cudaMemset(D_PW_ecoul_TOT, 0, sizeof(double));

    const int thread_num = 128;
    dim3 block_sz = dim3(thread_num);

    int tile_num_per_block = thread_num / 32;

    int x_blocks_num = (n_patoms + 31) / 32;           // solute
    int y_blocks_num = ((3 * n_waters) + 31) / 32;     // water

    int total_tiles = x_blocks_num * y_blocks_num;     // full rectangle

    int grid_sz = (total_tiles + tile_num_per_block - 1) / tile_num_per_block;

    dim3 grid = dim3(grid_sz);

    calc_pw<<<grid, block_sz>>>(
        n_patoms,
        3 * n_waters,
        D_charges,
        D_ccharges,
        D_atypes,
        D_catypes,
        topo,
        D_excluded,
        ctx.d_LJ_matrix,
        D_patoms,
        X,
        DV_X,
        n_atoms_solute,
        true,
        D_PW_evdw_TOT,
        D_PW_ecoul_TOT);

    cudaDeviceSynchronize();
    size_t mem_size_all = (n_atoms_solute + 3 * n_waters) * sizeof(dvel_t);
    cudaMemcpy(dvelocities, DV_X, mem_size_all, cudaMemcpyDeviceToHost);

    double PW_evdw_TOT, PW_ecoul_TOT;
    cudaMemcpy(&PW_evdw_TOT, D_PW_evdw_TOT, sizeof(double), cudaMemcpyDeviceToHost);
    cudaMemcpy(&PW_ecoul_TOT, D_PW_ecoul_TOT, sizeof(double), cudaMemcpyDeviceToHost);

    E_nonbond_pw.Uvdw += PW_evdw_TOT;
    E_nonbond_pw.Ucoul += PW_ecoul_TOT;
}

void init_nonbonded_pw_force_kernel_data() {
    using namespace CudaNonbondedPWForce;
    if (!is_initialized) {
        check_cudaMalloc((void**)&D_PW_evdw_TOT, sizeof(double));
        check_cudaMalloc((void**)&D_PW_ecoul_TOT, sizeof(double));

        is_initialized = true;
    }
}

void cleanup_nonbonded_pw_force() {
    using namespace CudaNonbondedPWForce;
    if (is_initialized) {
        cudaFree(D_PW_evdw_TOT);
        cudaFree(D_PW_ecoul_TOT);
        is_initialized = false;
    }
}
