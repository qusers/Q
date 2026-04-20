#pragma once

#include <cuda_runtime.h>

#include <array>
#include <map>
#include <vector>

#include "common/include/md_types.h"

#include "cuda_utility.cuh"

class CudaContext {
   public:
    /*
    Used in cuda_shake_constraints.cu
    */
    int* d_mol_n_shakes = nullptr;
    shake_bond_t* d_shake_bonds = nullptr;
    double* d_winv = nullptr;
    coord_t* d_xcoords = nullptr;

    /*
    Used in cuda_nonbonded_qq_force.cu
    */
    q_elscale_t* d_q_elscales = nullptr;
    atype_t* d_q_atypes = nullptr;
    E_nonbonded_t* d_EQ_nonbond_qq = nullptr;
    double* d_lambdas = nullptr;

    /*
    Used in cuda_polx_water_force.cu
    */
    shell_t* d_wshells = nullptr;
    /*
    Used in cuda_pshell_force.cu
    */
    bool* d_shell = nullptr;

    /*
    Used in cuda_restrang_force.cu
    */
    E_restraint_t* d_EQ_restraint = nullptr;

    /*
    Other helper arrays
    */
    ccharge_t* d_charge_table_all;  // Device copy of h_charge_table_all
    double* d_charge_pair_products = nullptr;
    catype_t* d_catype_table_all;   // Device copy of h_catype_table_all
    vdw_pair_param_t* d_catype_pair_params = nullptr;
    ccharge_t* d_unified_ccharges = nullptr;
    catype_t* d_unified_catypes = nullptr;
    int3* d_ngbrs_14 = nullptr;
    std::map<std::array<double, 4>, int> catype_to_type_host;
    int n_charge_types = 0;
    int zero_charge_type = -1;
    int n_catype_types = 0;
    int zero_catype_type = -1;
    int* d_p_charge_types;
    int* d_w_charge_types;
    int* d_q_charge_types;  // [0, lambdas * n_qatoms) is the normal q_charge_type, [lambdas * n_qatoms, ... ) is the lambda-scaled q_charge_type]
    
    int* d_p_catype_types;
    int* d_w_catype_types;
    int* d_q_catype_types;  // [0, lambdas * n_qatoms) is the normal q_catype_type, [lambdas * n_qatoms, ... ) is the lambda-scaled q_catype_type]

    int *d_p_atoms_list;
    int *d_w_atoms_list;
    int *d_q_atoms_list;
    std::vector<int> h_p_atoms_list;
    std::vector<int> h_w_atoms_list;
    std::vector<int> h_q_atoms_list;

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

    std::array<double, 4> get_catype_key(const catype_t& catype) {
        return {
            catype.aii_normal,
            catype.bii_normal,
            catype.aii_1_4,
            catype.bii_1_4};
    }

   private:
    CudaContext() = default;

   void free();

    ~CudaContext() { free(); }
    CudaContext(const CudaContext&) = delete;
    CudaContext& operator=(const CudaContext&) = delete;

    void initialize_charge_tables_host();
    void initialize_catype_tables_host();
    void initialize_atom_lists_host();
    void initialize_ngbrs14_host();
};
template <typename T>
void CudaContext::sync_array_to_device(T* dst, const T* src, int count) {
    check_cuda(cudaMemcpy(dst, src, count * sizeof(T), cudaMemcpyHostToDevice));
}
template <typename T>
void CudaContext::sync_array_to_host(T* dst, const T* src, int count) {
    check_cuda(cudaMemcpy(dst, src, count * sizeof(T), cudaMemcpyDeviceToHost));
}
