#pragma once

#include <cuda_runtime.h>

#include "system.h"
#include "utils.h"

class CudaContext {
   public:
    /*
    Common data
    */
    coord_t* d_coords = nullptr;
    dvel_t* d_dvelocities = nullptr;
    vel_t* d_velocities = nullptr;

    /*
    Used in CudaAngleForce.cu
    */
    angle_t* d_angles = nullptr;
    cangle_t* d_cangles = nullptr;

    /*
    Used in CudaBondForce.cu
    */
    bond_t* d_bonds = nullptr;
    cbond_t* d_cbonds = nullptr;

    /*
    Used in CudaImproper2Force.cu
    */
    improper_t* d_impropers = nullptr;
    cimproper_t* d_cimpropers = nullptr;

    /*
    Used in CudaLeapfrog.cu
    */
    atype_t* d_atypes = nullptr;
    catype_t* d_catypes = nullptr;

    /*
    Used in CudaShakeConstraints.cu
    */
    int* d_mol_n_shakes = nullptr;
    shake_bond_t* d_shake_bonds = nullptr;
    double* d_winv = nullptr;
    coord_t* d_xcoords = nullptr;

    /*
    Used in CudaNonbondedQQForce.cu
    */
    q_atom_t* d_q_atoms = nullptr;
    q_charge_t* d_q_charges = nullptr;
    int* d_LJ_matrix = nullptr;
    bool* d_excluded = nullptr;
    q_elscale_t* d_q_elscales = nullptr;
    q_catype_t* d_q_catypes = nullptr;
    q_atype_t* d_q_atypes = nullptr;
    E_nonbonded_t* d_EQ_nonbond_qq = nullptr;
    double* d_lambdas = nullptr;

    /*
    Used in CudaPolxWaterForce.cu
    */
    shell_t* d_wshells = nullptr;
    /*
    Used in CudaPshellForce.cu
    */
    bool* d_shell = nullptr;
    coord_t* d_coords_top = nullptr;

    /*
    Used in CudaRestrangForce.cu
    */
    restrang_t* d_restrangs = nullptr;
    E_restraint_t* d_EQ_restraint = nullptr;
    restrdis_t* d_restrdists = nullptr;

    /*
    Used in CudaRestrposForce.cu
    */
    restrpos_t* d_restrspos = nullptr;

    /*
    Used in CudaRestrseqForce.cu
    */
    restrseq_t* d_restrseqs = nullptr;
    bool* d_heavy = nullptr;

    /*
    Used in CudaRestrwallForce.cu
    */
    restrwall_t* d_restrwalls = nullptr;

    /*
    Used in CudaTorsionForce.cu
    */
    torsion_t* d_torsions = nullptr;
    ctorsion_t* d_ctorsions = nullptr;

    /*
    Used in CudaNonbondedPPForce.cu
    */
    ccharge_t* d_ccharges;
    charge_t* d_charges;
    p_atom_t* d_p_atoms;

    static CudaContext& instance() {
        static CudaContext ctx;
        return ctx;
    }
    void init();

    template <typename T>
    void sync_array_to_device(T* dst, const T* src, int count);

    template <typename T>
    void sync_array_to_host(T* dst, const T* src, int count);

    void sync_all_to_device();
    void sync_all_to_host();
    void reset_energies();

   private:
    CudaContext() = default;

    void free();

    ~CudaContext() { free(); }
    CudaContext(const CudaContext&) = delete;
    CudaContext& operator=(const CudaContext&) = delete;
};
template <typename T>
void CudaContext::sync_array_to_device(T* dst, const T* src, int count) {
    check_cuda(cudaMemcpy(dst, src, count * sizeof(T), cudaMemcpyHostToDevice));
}
template <typename T>
void CudaContext::sync_array_to_host(T* dst, const T* src, int count) {
    check_cuda(cudaMemcpy(dst, src, count * sizeof(T), cudaMemcpyDeviceToHost));
}