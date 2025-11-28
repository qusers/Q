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

    check_cudaMalloc((void**)&d_atypes, sizeof(atype_t) * n_atypes);
    check_cudaMalloc((void**)&d_catypes, sizeof(catype_t) * n_catypes);

    check_cudaMalloc((void**)&d_mol_shake_offset, sizeof(int) * n_molecules);  // calculation data, not initialized in the beginning.

    check_cudaMalloc((void**)&d_q_atoms, sizeof(q_atom_t) * n_qatoms);
    check_cudaMalloc((void**)&d_q_charges, sizeof(q_charge_t) * n_qatoms * n_lambdas);
    check_cudaMalloc((void**)&d_LJ_matrix, sizeof(int) * n_atoms_solute * n_atoms_solute);
    check_cudaMalloc((void**)&d_excluded, sizeof(bool) * n_atoms);
    check_cudaMalloc((void**)&d_q_elscales, sizeof(q_elscale_t) * n_qelscales);
    check_cudaMalloc((void**)&d_q_catypes, sizeof(q_catype_t) * n_qcatypes);
    check_cudaMalloc((void**)&d_q_atypes, sizeof(q_atype_t) * n_qatoms * n_lambdas);
    check_cudaMalloc((void**)&d_EQ_nonbond_qq, sizeof(E_nonbonded_t) * n_lambdas);
    check_cudaMalloc((void**)&d_lambdas, sizeof(double) * n_lambdas);

    check_cudaMalloc((void**)&d_wshells, n_shells * sizeof(shell_t));
    check_cudaMalloc((void**)&d_list_sh, n_max_inshell * n_shells * sizeof(int));  // calculation data, not initialized in the beginning.
    check_cudaMalloc((void**)&d_theta, n_waters * sizeof(double));                 // calculation data, not initialized in the beginning.
    check_cudaMalloc((void**)&d_theta0, n_waters * sizeof(double));                // calculation data, not initialized in the beginning.
    check_cudaMalloc((void**)&d_tdum, n_waters * sizeof(double));                  // calculation data, not initialized in the beginning.
    check_cudaMalloc((void**)&d_water_rank, n_waters * sizeof(int));               // calculation data, not initialized in the beginning.
    check_cudaMalloc((void**)&d_water_shell, n_waters * sizeof(int));              // calculation data, not initialized in the beginning.

    check_cudaMalloc((void**)&d_shell, sizeof(bool) * n_atoms);
    check_cudaMalloc((void**)&d_coords_top, sizeof(coord_t) * n_atoms);

    check_cudaMalloc((void**)&d_restrangs, sizeof(restrang_t) * n_restrangs);
    check_cudaMalloc((void**)&d_EQ_restraint, sizeof(E_restraint_t) * n_lambdas);
    check_cudaMalloc((void**)&d_restrdists, sizeof(restrdis_t) * n_restrdists);

    check_cudaMalloc((void**)&d_restrseqs, sizeof(restrseq_t) * n_restrseqs);
    check_cudaMalloc((void**)&d_heavy, sizeof(bool) * n_atoms);

    check_cudaMalloc((void**)&d_restrwalls, sizeof(restrwall_t) * n_restrwalls);

    check_cudaMalloc((void**)&d_torsions, sizeof(torsion_t) * n_torsions);
    check_cudaMalloc((void**)&d_ctorsions, sizeof(ctorsion_t) * n_ctorsions);

    check_cudaMalloc((void**)&d_ccharges, sizeof(ccharge_t) * n_ccharges);
    check_cudaMalloc((void**)&d_charges, sizeof(charge_t) * n_charges);
    check_cudaMalloc((void**)&d_p_atoms, sizeof(p_atom_t) * n_patoms);

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

    sync_array_to_device<atype_t>(d_atypes, atypes, n_atypes);
    sync_array_to_device<catype_t>(d_catypes, catypes, n_catypes);

    sync_array_to_device<q_atom_t>(d_q_atoms, q_atoms, n_qatoms);
    sync_array_to_device<q_charge_t>(d_q_charges, q_charges, n_qatoms * n_lambdas);
    sync_array_to_device<int>(d_LJ_matrix, LJ_matrix, n_atoms_solute * n_atoms_solute);
    sync_array_to_device<bool>(d_excluded, excluded, n_atoms);
    sync_array_to_device<q_elscale_t>(d_q_elscales, q_elscales, n_qelscales);
    sync_array_to_device<q_catype_t>(d_q_catypes, q_catypes, n_qcatypes);
    sync_array_to_device<q_atype_t>(d_q_atypes, q_atypes, n_qatoms * n_lambdas);
    sync_array_to_device<E_nonbonded_t>(d_EQ_nonbond_qq, EQ_nonbond_qq, n_lambdas);
    sync_array_to_device<double>(d_lambdas, lambdas, n_lambdas);
    sync_array_to_device<shell_t>(d_wshells, wshells, n_shells);

    sync_array_to_device<bool>(d_shell, shell, n_atoms);
    sync_array_to_device<coord_t>(d_coords_top, coords_top, n_atoms);

    sync_array_to_device<restrang_t>(d_restrangs, restrangs, n_restrangs);
    sync_array_to_device<E_restraint_t>(d_EQ_restraint, EQ_restraint, n_lambdas);
    sync_array_to_device<restrdis_t>(d_restrdists, restrdists, n_restrdists);

    sync_array_to_device<restrseq_t>(d_restrseqs, restrseqs, n_restrseqs);
    sync_array_to_device<bool>(d_heavy, heavy, n_atoms);
    sync_array_to_device<restrwall_t>(d_restrwalls, restrwalls, n_restrwalls);

    sync_array_to_device<torsion_t>(d_torsions, torsions, n_torsions);
    sync_array_to_device<ctorsion_t>(d_ctorsions, ctorsions, n_ctorsions);

    sync_array_to_device<ccharge_t>(d_ccharges, ccharges, n_ccharges);
    sync_array_to_device<charge_t>(d_charges, charges, n_charges);
    sync_array_to_device<p_atom_t>(d_p_atoms, p_atoms, n_patoms);
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

    sync_array_to_host<atype_t>(atypes, d_atypes, n_atypes);
    sync_array_to_host<catype_t>(catypes, d_catypes, n_catypes);

    sync_array_to_host<q_atom_t>(q_atoms, d_q_atoms, n_qatoms);
    sync_array_to_host<q_charge_t>(q_charges, d_q_charges, n_qatoms * n_lambdas);
    sync_array_to_host<int>(LJ_matrix, d_LJ_matrix, n_atoms_solute * n_atoms_solute);
    sync_array_to_host<bool>(excluded, d_excluded, n_atoms);
    sync_array_to_host<q_elscale_t>(q_elscales, d_q_elscales, n_qelscales);
    sync_array_to_host<q_catype_t>(q_catypes, d_q_catypes, n_qcatypes);
    sync_array_to_host<q_atype_t>(q_atypes, d_q_atypes, n_qatoms * n_lambdas);
    sync_array_to_host<E_nonbonded_t>(EQ_nonbond_qq, d_EQ_nonbond_qq, n_lambdas);
    sync_array_to_host<double>(lambdas, d_lambdas, n_lambdas);
    sync_array_to_host<shell_t>(wshells, d_wshells, n_shells);
    sync_array_to_host<bool>(shell, d_shell, n_atoms);
    sync_array_to_host<coord_t>(coords, d_coords_top, n_atoms);

    sync_array_to_host<restrang_t>(restrangs, d_restrangs, n_restrangs);
    sync_array_to_host<E_restraint_t>(EQ_restraint, d_EQ_restraint, n_lambdas);
    sync_array_to_host<restrdis_t>(restrdists, d_restrdists, n_restrdists);

    sync_array_to_host<restrseq_t>(restrseqs, d_restrseqs, n_restrseqs);
    sync_array_to_host<bool>(heavy, d_heavy, n_atoms);
    sync_array_to_host<restrwall_t>(restrwalls, d_restrwalls, n_restrwalls);

    sync_array_to_host<torsion_t>(torsions, d_torsions, n_torsions);
    sync_array_to_host<ctorsion_t>(ctorsions, d_ctorsions, n_ctorsions);

    sync_array_to_host<ccharge_t>(ccharges, d_ccharges, n_ccharges);
    sync_array_to_host<charge_t>(charges, d_charges, n_charges);
    sync_array_to_host<p_atom_t>(p_atoms, d_p_atoms, n_patoms);
}

void CudaContext::free() {
}
