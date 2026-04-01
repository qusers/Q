#include "cuda_context.cuh"
#include "context.h"

void CudaContext::init() {
    auto& host = Context::instance();

    check_cudaMalloc((void**)&d_coords, sizeof(coord_t) * host.n_atoms);
    check_cudaMalloc((void**)&d_dvelocities, sizeof(dvel_t) * host.n_atoms);
    check_cudaMalloc((void**)&d_velocities, sizeof(vel_t) * host.n_atoms);
    check_cudaMalloc((void**)&d_angles, sizeof(angle_t) * host.n_angles);
    check_cudaMalloc((void**)&d_cangles, sizeof(cangle_t) * host.n_cangles);
    check_cudaMalloc((void**)&d_bonds, sizeof(bond_t) * host.n_bonds);
    check_cudaMalloc((void**)&d_cbonds, sizeof(cbond_t) * host.n_cbonds);
    check_cudaMalloc((void**)&d_impropers, sizeof(improper_t) * host.n_impropers);
    check_cudaMalloc((void**)&d_cimpropers, sizeof(cimproper_t) * host.n_cimpropers);

    check_cudaMalloc((void**)&d_mol_n_shakes, sizeof(int) * host.n_molecules);
    check_cudaMalloc((void**)&d_shake_bonds, sizeof(shake_bond_t) * host.n_shake_constraints);
    check_cudaMalloc((void**)&d_winv, sizeof(double) * host.n_atoms);
    check_cudaMalloc((void**)&d_xcoords, sizeof(coord_t) * host.n_atoms);

    check_cudaMalloc((void**)&d_atypes, sizeof(atype_t) * host.n_atypes);
    check_cudaMalloc((void**)&d_catypes, sizeof(catype_t) * host.n_catypes);

    check_cudaMalloc((void**)&d_q_atoms, sizeof(q_atom_t) * host.n_qatoms);
    check_cudaMalloc((void**)&d_q_charges, sizeof(q_charge_t) * host.n_qatoms * host.n_lambdas);
    check_cudaMalloc((void**)&d_LJ_matrix, sizeof(int) * host.n_atoms_solute * host.n_atoms_solute);
    check_cudaMalloc((void**)&d_excluded, sizeof(bool) * host.n_atoms);
    check_cudaMalloc((void**)&d_q_elscales, sizeof(q_elscale_t) * host.n_qelscales);
    check_cudaMalloc((void**)&d_q_catypes, sizeof(q_catype_t) * host.n_qcatypes);
    check_cudaMalloc((void**)&d_q_atypes, sizeof(q_atype_t) * host.n_qatoms * host.n_lambdas);
    check_cudaMalloc((void**)&d_EQ_nonbond_qq, sizeof(E_nonbonded_t) * host.n_lambdas);
    check_cudaMalloc((void**)&d_lambdas, sizeof(double) * host.n_lambdas);

    check_cudaMalloc((void**)&d_wshells, host.n_shells * sizeof(shell_t));
    check_cudaMalloc((void**)&d_shell, sizeof(bool) * host.n_atoms);
    check_cudaMalloc((void**)&d_coords_top, sizeof(coord_t) * host.n_atoms);

    check_cudaMalloc((void**)&d_restrangs, sizeof(restrang_t) * host.n_restrangs);
    check_cudaMalloc((void**)&d_EQ_restraint, sizeof(E_restraint_t) * host.n_lambdas);
    check_cudaMalloc((void**)&d_restrdists, sizeof(restrdis_t) * host.n_restrdists);

    check_cudaMalloc((void**)&d_restrspos, sizeof(restrpos_t) * host.n_restrspos);

    check_cudaMalloc((void**)&d_restrseqs, sizeof(restrseq_t) * host.n_restrseqs);
    check_cudaMalloc((void**)&d_heavy, sizeof(bool) * host.n_atoms);

    check_cudaMalloc((void**)&d_restrwalls, sizeof(restrwall_t) * host.n_restrwalls);

    check_cudaMalloc((void**)&d_torsions, sizeof(torsion_t) * host.n_torsions);
    check_cudaMalloc((void**)&d_ctorsions, sizeof(ctorsion_t) * host.n_ctorsions);

    check_cudaMalloc((void**)&d_ccharges, sizeof(ccharge_t) * host.n_ccharges);
    check_cudaMalloc((void**)&d_charges, sizeof(charge_t) * host.n_charges);
    check_cudaMalloc((void**)&d_p_atoms, sizeof(p_atom_t) * host.n_patoms);

    initialize_atom_lists_host();
    initialize_charge_tables_host();
    initialize_catype_tables_host();

    sync_all_to_device();
}

void CudaContext::sync_all_to_device() {
    auto& host = Context::instance();

    sync_array_to_device<coord_t>(d_coords, host.coords.data(), host.n_atoms);
    sync_array_to_device<dvel_t>(d_dvelocities, host.dvelocities.data(), host.n_atoms);
    sync_array_to_device<vel_t>(d_velocities, host.velocities.data(), host.n_atoms);
    sync_array_to_device<angle_t>(d_angles, host.angles.data(), host.n_angles);
    sync_array_to_device<cangle_t>(d_cangles, host.cangles.data(), host.n_cangles);
    sync_array_to_device<bond_t>(d_bonds, host.bonds.data(), host.n_bonds);
    sync_array_to_device<cbond_t>(d_cbonds, host.cbonds.data(), host.n_cbonds);
    sync_array_to_device<improper_t>(d_impropers, host.impropers.data(), host.n_impropers);
    sync_array_to_device<cimproper_t>(d_cimpropers, host.cimpropers.data(), host.n_cimpropers);

    sync_array_to_device<int>(d_mol_n_shakes, host.mol_n_shakes.data(), host.n_molecules);
    sync_array_to_device<shake_bond_t>(d_shake_bonds, host.shake_bonds.data(), host.n_shake_constraints);
    sync_array_to_device<double>(d_winv, host.winv.data(), host.n_atoms);
    sync_array_to_device<coord_t>(d_xcoords, host.xcoords.data(), host.n_atoms);

    sync_array_to_device<atype_t>(d_atypes, host.atypes.data(), host.n_atypes);
    sync_array_to_device<catype_t>(d_catypes, host.catypes.data(), host.n_catypes);

    sync_array_to_device<q_atom_t>(d_q_atoms, host.q_atoms.data(), host.n_qatoms);
    sync_array_to_device<q_charge_t>(d_q_charges, host.q_charges.data(), host.n_qatoms * host.n_lambdas);
    sync_array_to_device<int>(d_LJ_matrix, host.LJ_matrix.data(), host.n_atoms_solute * host.n_atoms_solute);
    sync_array_to_device<bool>(d_excluded, host.excluded.get(), host.n_atoms);
    sync_array_to_device<q_elscale_t>(d_q_elscales, host.q_elscales.data(), host.n_qelscales);
    sync_array_to_device<q_catype_t>(d_q_catypes, host.q_catypes.data(), host.n_qcatypes);
    sync_array_to_device<q_atype_t>(d_q_atypes, host.q_atypes.data(), host.n_qatoms * host.n_lambdas);
    sync_array_to_device<E_nonbonded_t>(d_EQ_nonbond_qq, host.EQ_nonbond_qq.data(), host.n_lambdas);
    sync_array_to_device<double>(d_lambdas, host.lambdas.data(), host.n_lambdas);
    sync_array_to_device<shell_t>(d_wshells, host.wshells.data(), host.n_shells);

    sync_array_to_device<bool>(d_shell, host.shell.get(), host.n_atoms);
    sync_array_to_device<coord_t>(d_coords_top, host.coords_top.data(), host.n_atoms);

    sync_array_to_device<restrang_t>(d_restrangs, host.restrangs.data(), host.n_restrangs);
    sync_array_to_device<E_restraint_t>(d_EQ_restraint, host.EQ_restraint.data(), host.n_lambdas);
    sync_array_to_device<restrdis_t>(d_restrdists, host.restrdists.data(), host.n_restrdists);

    sync_array_to_device<restrpos_t>(d_restrspos, host.restrspos.data(), host.n_restrspos);

    sync_array_to_device<restrseq_t>(d_restrseqs, host.restrseqs.data(), host.n_restrseqs);
    sync_array_to_device<bool>(d_heavy, host.heavy.get(), host.n_atoms);
    sync_array_to_device<restrwall_t>(d_restrwalls, host.restrwalls.data(), host.n_restrwalls);

    sync_array_to_device<torsion_t>(d_torsions, host.torsions.data(), host.n_torsions);
    sync_array_to_device<ctorsion_t>(d_ctorsions, host.ctorsions.data(), host.n_ctorsions);

    sync_array_to_device<ccharge_t>(d_ccharges, host.ccharges.data(), host.n_ccharges);
    sync_array_to_device<charge_t>(d_charges, host.charges.data(), host.n_charges);
    sync_array_to_device<p_atom_t>(d_p_atoms, host.p_atoms.data(), host.n_patoms);
}

void CudaContext::sync_all_to_host() {
    auto& host = Context::instance();

    sync_array_to_host<coord_t>(host.coords.data(), d_coords, host.n_atoms);
    sync_array_to_host<dvel_t>(host.dvelocities.data(), d_dvelocities, host.n_atoms);
    sync_array_to_host<vel_t>(host.velocities.data(), d_velocities, host.n_atoms);
    sync_array_to_host<angle_t>(host.angles.data(), d_angles, host.n_angles);
    sync_array_to_host<cangle_t>(host.cangles.data(), d_cangles, host.n_cangles);
    sync_array_to_host<bond_t>(host.bonds.data(), d_bonds, host.n_bonds);
    sync_array_to_host<cbond_t>(host.cbonds.data(), d_cbonds, host.n_cbonds);
    sync_array_to_host<improper_t>(host.impropers.data(), d_impropers, host.n_impropers);
    sync_array_to_host<cimproper_t>(host.cimpropers.data(), d_cimpropers, host.n_cimpropers);

    sync_array_to_host<int>(host.mol_n_shakes.data(), d_mol_n_shakes, host.n_molecules);
    sync_array_to_host<shake_bond_t>(host.shake_bonds.data(), d_shake_bonds, host.n_shake_constraints);
    sync_array_to_host<double>(host.winv.data(), d_winv, host.n_atoms);
    sync_array_to_host<coord_t>(host.xcoords.data(), d_xcoords, host.n_atoms);

    sync_array_to_host<atype_t>(host.atypes.data(), d_atypes, host.n_atypes);
    sync_array_to_host<catype_t>(host.catypes.data(), d_catypes, host.n_catypes);

    sync_array_to_host<q_atom_t>(host.q_atoms.data(), d_q_atoms, host.n_qatoms);
    sync_array_to_host<q_charge_t>(host.q_charges.data(), d_q_charges, host.n_qatoms * host.n_lambdas);
    sync_array_to_host<int>(host.LJ_matrix.data(), d_LJ_matrix, host.n_atoms_solute * host.n_atoms_solute);
    sync_array_to_host<bool>(host.excluded.get(), d_excluded, host.n_atoms);
    sync_array_to_host<q_elscale_t>(host.q_elscales.data(), d_q_elscales, host.n_qelscales);
    sync_array_to_host<q_catype_t>(host.q_catypes.data(), d_q_catypes, host.n_qcatypes);
    sync_array_to_host<q_atype_t>(host.q_atypes.data(), d_q_atypes, host.n_qatoms * host.n_lambdas);
    sync_array_to_host<E_nonbonded_t>(host.EQ_nonbond_qq.data(), d_EQ_nonbond_qq, host.n_lambdas);
    sync_array_to_host<double>(host.lambdas.data(), d_lambdas, host.n_lambdas);
    sync_array_to_host<shell_t>(host.wshells.data(), d_wshells, host.n_shells);
    sync_array_to_host<bool>(host.shell.get(), d_shell, host.n_atoms);
    sync_array_to_host<coord_t>(host.coords_top.data(), d_coords_top, host.n_atoms);

    sync_array_to_host<restrang_t>(host.restrangs.data(), d_restrangs, host.n_restrangs);
    sync_array_to_host<E_restraint_t>(host.EQ_restraint.data(), d_EQ_restraint, host.n_lambdas);
    sync_array_to_host<restrdis_t>(host.restrdists.data(), d_restrdists, host.n_restrdists);

    sync_array_to_host<restrpos_t>(host.restrspos.data(), d_restrspos, host.n_restrspos);

    sync_array_to_host<restrseq_t>(host.restrseqs.data(), d_restrseqs, host.n_restrseqs);
    sync_array_to_host<bool>(host.heavy.get(), d_heavy, host.n_atoms);
    sync_array_to_host<restrwall_t>(host.restrwalls.data(), d_restrwalls, host.n_restrwalls);

    sync_array_to_host<torsion_t>(host.torsions.data(), d_torsions, host.n_torsions);
    sync_array_to_host<ctorsion_t>(host.ctorsions.data(), d_ctorsions, host.n_ctorsions);

    sync_array_to_host<ccharge_t>(host.ccharges.data(), d_ccharges, host.n_ccharges);
    sync_array_to_host<charge_t>(host.charges.data(), d_charges, host.n_charges);
    sync_array_to_host<p_atom_t>(host.p_atoms.data(), d_p_atoms, host.n_patoms);
}

void CudaContext::free() {
    cudaFree(d_coords);
    cudaFree(d_dvelocities);
    cudaFree(d_velocities);
    cudaFree(d_angles);
    cudaFree(d_cangles);
    cudaFree(d_bonds);
    cudaFree(d_cbonds);
    cudaFree(d_impropers);
    cudaFree(d_cimpropers);

    cudaFree(d_atypes);
    cudaFree(d_catypes);

    cudaFree(d_mol_n_shakes);
    cudaFree(d_shake_bonds);
    cudaFree(d_winv);
    cudaFree(d_xcoords);

    cudaFree(d_q_atoms);
    cudaFree(d_q_charges);
    cudaFree(d_LJ_matrix);
    cudaFree(d_excluded);
    cudaFree(d_q_elscales);
    cudaFree(d_q_catypes);
    cudaFree(d_q_atypes);
    cudaFree(d_EQ_nonbond_qq);
    cudaFree(d_lambdas);

    cudaFree(d_wshells);
    cudaFree(d_shell);
    cudaFree(d_coords_top);

    cudaFree(d_restrangs);
    cudaFree(d_EQ_restraint);
    cudaFree(d_restrdists);

    cudaFree(d_restrspos);

    cudaFree(d_restrseqs);
    cudaFree(d_heavy);

    cudaFree(d_restrwalls);

    cudaFree(d_torsions);
    cudaFree(d_ctorsions);

    cudaFree(d_ccharges);
    cudaFree(d_charges);
    cudaFree(d_p_atoms);
    cudaFree(d_charge_table_all);
    cudaFree(d_catype_table_all);
    cudaFree(d_p_charge_types);
    cudaFree(d_w_charge_types);
    cudaFree(d_q_charge_types);
    cudaFree(d_p_catype_types);
    cudaFree(d_w_catype_types);
    cudaFree(d_q_catype_types);
    cudaFree(d_p_atoms_list);
    cudaFree(d_w_atoms_list);
    cudaFree(d_q_atoms_list);

    d_charge_table_all = nullptr;
    d_catype_table_all = nullptr;
    d_p_charge_types = nullptr;
    d_w_charge_types = nullptr;
    d_q_charge_types = nullptr;
    d_p_catype_types = nullptr;
    d_w_catype_types = nullptr;
    d_q_catype_types = nullptr;
    d_p_atoms_list = nullptr;
    d_w_atoms_list = nullptr;
    d_q_atoms_list = nullptr;
}

void CudaContext::reset_energies() {
    auto& host = Context::instance();
    cudaMemset(d_dvelocities, 0, sizeof(dvel_t) * host.n_atoms);
    cudaMemset(d_EQ_nonbond_qq, 0, sizeof(E_nonbonded_t) * host.n_lambdas);
    cudaMemset(d_EQ_restraint, 0, sizeof(E_restraint_t) * host.n_lambdas);
}

void CudaContext::initialize_atom_lists_host() {
    auto& host = Context::instance();

    h_p_atoms_list.clear();
    for (int i = 0; i < host.n_patoms; i++) {
        int id = host.p_atoms[i].a - 1;
        if (!host.excluded[id]) {
            h_p_atoms_list.push_back(id);
        }
    }

    h_w_atoms_list.clear();
    for (int i = host.n_atoms_solute; i < host.n_atoms; i++) {
        int id = i;
        if (!host.excluded[id]) {
            h_w_atoms_list.push_back(id);
        }
    }
    printf("Number of water atoms: %d, number of all water atoms: %d\n", (int)h_w_atoms_list.size(), host.n_atoms - host.n_atoms_solute);

    h_q_atoms_list.clear();
    for (int i = 0; i < host.n_qatoms; i++) {
        int id = host.q_atoms[i].a - 1;
        if (!host.excluded[id]) {
            h_q_atoms_list.push_back(id);
        }
    }

    check_cudaMalloc((void**)&d_p_atoms_list, sizeof(int) * h_p_atoms_list.size());
    cudaMemcpy(d_p_atoms_list, h_p_atoms_list.data(), sizeof(int) * h_p_atoms_list.size(), cudaMemcpyHostToDevice);
    check_cudaMalloc((void**)&d_w_atoms_list, sizeof(int) * h_w_atoms_list.size());
    cudaMemcpy(d_w_atoms_list, h_w_atoms_list.data(), sizeof(int) * h_w_atoms_list.size(), cudaMemcpyHostToDevice);
    check_cudaMalloc((void**)&d_q_atoms_list, sizeof(int) * h_q_atoms_list.size());
    cudaMemcpy(d_q_atoms_list, h_q_atoms_list.data(), sizeof(int) * h_q_atoms_list.size(), cudaMemcpyHostToDevice);
}

void CudaContext::initialize_charge_tables_host() {
    auto& host = Context::instance();

    std::map<double, int> charge_to_type_host;
    std::vector<ccharge_t> h_charge_table_all;  // h_charge_table_all[charge type] = (code, charge)

    auto add_charge = [&](double charge) {
        if (charge_to_type_host.count(charge) == 0) {
            int sz = h_charge_table_all.size();
            ccharge_t new_ccharge;
            new_ccharge.code = sz;
            new_ccharge.charge = charge;
            h_charge_table_all.push_back(new_ccharge);
            charge_to_type_host[charge] = sz;
        }
    };

    for (int i = 0; i < host.n_ccharges; i++) {
        double charge = host.ccharges[i].charge;
        add_charge(charge);
    }
    for (int state = 0; state < host.n_lambdas; state++) {
        for (int i = 0; i < host.n_qatoms; i++) {
            double charge = host.q_charges[i + host.n_qatoms * state].q;
            add_charge(charge);
            charge *= host.lambdas[state];
            add_charge(charge);
        }
    }

    std::vector<int> p_charge_types(h_p_atoms_list.size());
    for (int i = 0; i < static_cast<int>(h_p_atoms_list.size()); i++) {
        int id = h_p_atoms_list[i];
        double charge = host.ccharges[host.charges[id].code - 1].charge;
        p_charge_types[i] = charge_to_type_host[charge];
    }

    int h_q_atoms_size = h_q_atoms_list.size();
    std::vector<int> q_charge_types(h_q_atoms_size * host.n_lambdas);
    std::map<int, int> q_idx;
    for (int i = 0; i < host.n_qatoms; i++) {
        int id = host.q_atoms[i].a - 1;
        if (!host.excluded[id]) {
            q_idx[id] = i;
        }
    }

    for (int state = 0; state < host.n_lambdas; state++) {
        for (int i = 0; i < h_q_atoms_size; i++) {
            int id = h_q_atoms_list[i];
            double charge = host.q_charges[q_idx[id] + host.n_qatoms * state].q;
            q_charge_types[state * h_q_atoms_size + i] = charge_to_type_host[charge];
        }
    }

    int h_w_atoms_size = h_w_atoms_list.size();
    std::vector<int> w_charge_types(h_w_atoms_size);
    for (int i = 0; i < h_w_atoms_size; i++) {
        int id = h_w_atoms_list[i];
        double charge = host.ccharges[host.charges[id].code - 1].charge;
        w_charge_types[i] = charge_to_type_host[charge];
    }

    check_cudaMalloc((void**)&d_charge_table_all, sizeof(ccharge_t) * h_charge_table_all.size());
    cudaMemcpy(d_charge_table_all, h_charge_table_all.data(), sizeof(ccharge_t) * h_charge_table_all.size(), cudaMemcpyHostToDevice);
    check_cudaMalloc((void**)&d_p_charge_types, sizeof(int) * p_charge_types.size());
    cudaMemcpy(d_p_charge_types, p_charge_types.data(), sizeof(int) * p_charge_types.size(), cudaMemcpyHostToDevice);
    check_cudaMalloc((void**)&d_q_charge_types, sizeof(int) * q_charge_types.size());
    cudaMemcpy(d_q_charge_types, q_charge_types.data(), sizeof(int) * q_charge_types.size(), cudaMemcpyHostToDevice);
    check_cudaMalloc((void**)&d_w_charge_types, sizeof(int) * w_charge_types.size());
    cudaMemcpy(d_w_charge_types, w_charge_types.data(), sizeof(int) * w_charge_types.size(), cudaMemcpyHostToDevice);
}

void CudaContext::initialize_catype_tables_host() {
    auto& host = Context::instance();

    std::vector<catype_t> h_catype_table_all;  // h_catype_table_all[catype code] = catype_t
    catype_to_type_host.clear();

    auto add_catype = [&](catype_t catype) {
        auto key = get_catype_key(catype);
        if (catype_to_type_host.count(key) == 0) {
            int sz = h_catype_table_all.size();
            catype.code = sz;
            h_catype_table_all.push_back(catype);
            catype_to_type_host[key] = sz;
        }
    };

    for (int i = 0; i < host.n_catypes; i++) {
        add_catype(host.catypes[i]);
    }

    for (int state = 0; state < host.n_lambdas; state++) {
        for (int i = 0; i < host.n_qatoms; i++) {
            const q_atype_t& qat = host.q_atypes[i + host.n_qatoms * state];
            const q_catype_t& catype = host.q_catypes[qat.code - 1];

            catype_t new_catype;
            new_catype.m = catype.m;
            new_catype.aii_normal = catype.Ai;
            new_catype.bii_normal = catype.Bi;
            new_catype.aii_polar = catype.Ci;
            new_catype.bii_polar = catype.ai;
            new_catype.aii_1_4 = catype.Ai_14;
            new_catype.bii_1_4 = catype.Bi_14;

            add_catype(new_catype);
        }
    }

    std::vector<int> p_catype_types(h_p_atoms_list.size());
    for (int i = 0; i < static_cast<int>(h_p_atoms_list.size()); i++) {
        int id = h_p_atoms_list[i];
        auto catype = host.catypes[host.atypes[id].code - 1];
        auto key = get_catype_key(catype);
        p_catype_types[i] = catype_to_type_host[key];
    }

    int h_q_atoms_size = h_q_atoms_list.size();
    std::vector<int> q_catype_types(h_q_atoms_size * host.n_lambdas);
    std::map<int, int> q_idx;
    for (int i = 0; i < host.n_qatoms; i++) {
        int id = host.q_atoms[i].a - 1;
        if (!host.excluded[id]) {
            q_idx[id] = i;
        }
    }


    for (int state = 0; state < host.n_lambdas; state++) {
        for (int i = 0; i < h_q_atoms_size; i++) {
            int id = h_q_atoms_list[i];
            const q_atype_t& qat = host.q_atypes[q_idx[id] + host.n_qatoms * state];
            const q_catype_t& catype = host.q_catypes[qat.code - 1];

            catype_t normal_catype;
            normal_catype.m = catype.m;
            normal_catype.aii_normal = catype.Ai;
            normal_catype.bii_normal = catype.Bi;
            normal_catype.aii_polar = catype.Ci;
            normal_catype.bii_polar = catype.ai;
            normal_catype.aii_1_4 = catype.Ai_14;
            normal_catype.bii_1_4 = catype.Bi_14;

            auto key_normal = get_catype_key(normal_catype);
            q_catype_types[state * h_q_atoms_size + i] = catype_to_type_host[key_normal];
        }
    }

    int h_w_atoms_size = h_w_atoms_list.size();
    std::vector<int> w_catype_types(h_w_atoms_size);
    for (int i = 0; i < h_w_atoms_size; i++) {
        int id = h_w_atoms_list[i];
        auto catype = host.catypes[host.atypes[id].code - 1];
        auto key = get_catype_key(catype);
        w_catype_types[i] = catype_to_type_host[key];
        // printf("Water atom %d: catype code %d mapped to %d, data: m=%f, aii_normal=%f, bii_normal=%f, aii_polar=%f, bii_polar=%f, aii_1_4=%f, bii_1_4=%f\n", id, atypes[id].code, w_catype_types[i], catype.m, catype.aii_normal, catype.bii_normal, catype.aii_polar, catype.bii_polar, catype.aii_1_4, catype.bii_1_4);
    }
    printf("Total water atom number: %d, w_catype_types size: %lu\n", h_w_atoms_size, w_catype_types.size());

    check_cudaMalloc((void**)&d_catype_table_all, sizeof(catype_t) * h_catype_table_all.size());
    cudaMemcpy(d_catype_table_all, h_catype_table_all.data(), sizeof(catype_t) * h_catype_table_all.size(), cudaMemcpyHostToDevice);
    check_cudaMalloc((void**)&d_p_catype_types, sizeof(int) * p_catype_types.size());
    cudaMemcpy(d_p_catype_types, p_catype_types.data(), sizeof(int) * p_catype_types.size(), cudaMemcpyHostToDevice);
    check_cudaMalloc((void**)&d_q_catype_types, sizeof(int) * q_catype_types.size());
    cudaMemcpy(d_q_catype_types, q_catype_types.data(), sizeof(int) * q_catype_types.size(), cudaMemcpyHostToDevice);
    check_cudaMalloc((void**)&d_w_catype_types, sizeof(int) * w_catype_types.size());
    cudaMemcpy(d_w_catype_types, w_catype_types.data(), sizeof(int) * w_catype_types.size(), cudaMemcpyHostToDevice);
}
