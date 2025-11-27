#include <iostream>

#include "cuda/include/CudaContext.cuh"

void CudaContext::init() {
    check_cudaMalloc((void**)&d_coords, sizeof(coord_t) * n_atoms);
    check_cudaMalloc((void**)&d_dvelocities, sizeof(dvel_t) * n_atoms);
    check_cudaMalloc((void**)&d_velocities, sizeof(vel_t) * n_atoms);
    check_cudaMalloc((void**)&d_angles, sizeof(angle_t) * n_angles);
    check_cudaMalloc((void**)&d_cangles, sizeof(cangle_t) * n_cangles);
    check_cudaMalloc((void**)&d_bonds, sizeof(bond_t) * n_bonds);
    check_cudaMalloc((void**)&d_cbonds, sizeof(cbond_t) * n_cbonds);
    check_cudaMalloc((void**)&d_impropers, sizeof(improper_t) * n_impropers);
    check_cudaMalloc((void**)&d_cimpropers, sizeof(cimproper_t) * n_cimpropers);

    check_cudaMalloc((void**)&d_mol_n_shakes, sizeof(int) * n_molecules);
    check_cudaMalloc((void**)&d_shake_bonds, sizeof(shake_bond_t) * n_shake_constraints);
    check_cudaMalloc((void**)&d_winv, sizeof(double) * n_atoms);
    check_cudaMalloc((void**)&d_xcoords, sizeof(coord_t) * n_atoms);

    check_cudaMalloc((void**)&d_mol_shake_offset, sizeof(int) * n_molecules);  // calculation data, not initialized in the beginning.

    sync_all_to_device();
}

void CudaContext::sync_all_to_device() {
    sync_array_to_device<coord_t>(d_coords, coords, n_atoms);
    sync_array_to_device<dvel_t>(d_dvelocities, dvelocities, n_atoms);
    sync_array_to_device<vel_t>(d_velocities, velocities, n_atoms);
    sync_array_to_device<angle_t>(d_angles, angles, n_angles);
    sync_array_to_device<cangle_t>(d_cangles, cangles, n_cangles);
    sync_array_to_device<bond_t>(d_bonds, bonds, n_bonds);
    sync_array_to_device<cbond_t>(d_cbonds, cbonds, n_cbonds);
    sync_array_to_device<improper_t>(d_impropers, impropers, n_impropers);
    sync_array_to_device<cimproper_t>(d_cimpropers, cimpropers, n_cimpropers);

    sync_array_to_device<int>(d_mol_n_shakes, mol_n_shakes, n_molecules);
    sync_array_to_device<shake_bond_t>(d_shake_bonds, shake_bonds, n_shake_constraints);
    sync_array_to_device<double>(d_winv, winv, n_atoms);
    sync_array_to_device<coord_t>(d_xcoords, xcoords, n_atoms);
}

void CudaContext::sync_all_to_host() {
    sync_array_to_host<coord_t>(coords, d_coords, n_atoms);
    sync_array_to_host<dvel_t>(dvelocities, d_dvelocities, n_atoms);
    sync_array_to_host<vel_t>(velocities, d_velocities, n_atoms);
    sync_array_to_host<angle_t>(angles, d_angles, n_angles);
    sync_array_to_host<cangle_t>(cangles, d_cangles, n_cangles);
    sync_array_to_host<bond_t>(bonds, d_bonds, n_bonds);
    sync_array_to_host<cbond_t>(cbonds, d_cbonds, n_cbonds);
    sync_array_to_host<improper_t>(impropers, d_impropers, n_impropers);
    sync_array_to_host<cimproper_t>(cimpropers, d_cimpropers, n_cimpropers);

    sync_array_to_host<int>(mol_n_shakes, d_mol_n_shakes, n_molecules);
    sync_array_to_host<shake_bond_t>(shake_bonds, d_shake_bonds, n_shake_constraints);
    sync_array_to_host<double>(winv, d_winv, n_atoms);
    sync_array_to_host<coord_t>(xcoords, d_xcoords, n_atoms);
}

void CudaContext::free() {
}
