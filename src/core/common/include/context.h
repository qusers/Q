#pragma once

#include <vector_types.h>

#include <array>
#include <memory>
#include <string>
#include <vector>

#include "common/include/md_types.h"
#include "common/include/nonbonded_14_mode.h"
#include "host_device_buffer.h"

class Context {
   public:
    static Context& instance() {
        static Context ctx;
        return ctx;
    }
    /* =============================================
     * == CONFIG
     * =============================================
     */

    std::string base_folder;
    bool run_gpu = false;

    int n_atoms = 0;         // the total number of atoms
    int n_atoms_solute = 0;  // the total number of solute number, in our system [0, n_atoms_solute) are solute, [n_atoms_solute, n_atoms) are water atoms
    int n_patoms = 0;
    int n_qatoms = 0;
    int n_waters = 0;
    int n_molecules = 0;
    double dt = 0.0;
    double tau_T = 0.0;
    md_t md;

    /*
    
    */
    std::unique_ptr<HostDeviceBuffer<coord_t>> coords;
    std::unique_ptr<HostDeviceBuffer<vel_t>> velocities;
    std::unique_ptr<HostDeviceBuffer<dvel_t>> dvelocities;


    /*
    */
    int n_angles = 0;
    int n_angles_solute = 0;
    int n_cangles = 0;
    std::unique_ptr<HostDeviceBuffer<angle_t>> angles;
    std::unique_ptr<HostDeviceBuffer<cangle_t>> cangles;

    int n_bonds = 0;
    int n_bonds_solute = 0;
    int n_cbonds = 0;
    std::unique_ptr<HostDeviceBuffer<bond_t>> bonds;
    std::unique_ptr<HostDeviceBuffer<cbond_t>> cbonds;
    int n_torsions = 0;
    int n_torsions_solute = 0;
    int n_ctorsions = 0;
    std::unique_ptr<HostDeviceBuffer<torsion_t>> torsions;
    std::unique_ptr<HostDeviceBuffer<ctorsion_t>> ctorsions;
    int n_impropers = 0;
    int n_impropers_solute = 0;
    int n_cimpropers = 0;
    std::unique_ptr<HostDeviceBuffer<improper_t>> impropers;
    std::unique_ptr<HostDeviceBuffer<cimproper_t>> cimpropers;

    int n_restrspos = 0;
    std::unique_ptr<HostDeviceBuffer<restrpos_t>> restrspos;

    int n_restrangs = 0;
    std::unique_ptr<HostDeviceBuffer<restrang_t>> restrangs;
    

    int n_restrdists = 0;
    std::unique_ptr<HostDeviceBuffer<restrdis_t>> restrdists;

    int n_restrseqs = 0;
    std::unique_ptr<HostDeviceBuffer<restrseq_t>> restrseqs;
    int n_restrwalls = 0;
    std::unique_ptr<HostDeviceBuffer<restrwall_t>> restrwalls;

    /*
    Atom Info
    */
    int n_charges = 0;
    int n_ccharges = 0;
    std::unique_ptr<HostDeviceBuffer<charge_t>> charges;
    std::unique_ptr<HostDeviceBuffer<ccharge_t>> ccharges;
    int n_atypes = 0;
    int n_catypes = 0;
    std::unique_ptr<HostDeviceBuffer<atype_t>> atypes;
    std::unique_ptr<HostDeviceBuffer<catype_t>> catypes;

    std::unique_ptr<HostDeviceBuffer<bool>> heavy;
    std::unique_ptr<HostDeviceBuffer<coord_t>> coords_init;

    std::unique_ptr<HostDeviceBuffer<bool>> excluded;

    std::unique_ptr<HostDeviceBuffer<double>> winv;

    std::unique_ptr<HostDeviceBuffer<bool>> shell;
    /*
    Pair
    */
    std::unique_ptr<HostDeviceBuffer<int>> LJ_matrix;

    /*
    Shake
    */
    int n_shake_constraints = 0;
    std::unique_ptr<HostDeviceBuffer<int>> mol_n_shakes;
    std::unique_ptr<HostDeviceBuffer<shake_bond_t>> shake_bonds;
    std::unique_ptr<HostDeviceBuffer<coord_t>> xcoords; // todo: It's just a temporary variables...
    /*
    Water
    */

    std::unique_ptr<HostDeviceBuffer<shell_t>> wshells;



    /*
    FEP
    */
    std::unique_ptr<HostDeviceBuffer<double>> lambdas; // Actually length is only 2..

    /*
    Energy
    */

    std::unique_ptr<HostDeviceBuffer<E_restraint_t>> EQ_restraint;

    /*
    */

    int n_ngbrs23 = 0;
    int n_ngbrs14 = 0;
    int n_excluded = 0;
    int n_cgrps_solute = 0;
    int n_cgrps_solvent = 0;
    int iuse_switch_atom = 0;

    std::vector<int> atom_to_qi;
    std::vector<ccharge_t> unified_ccharges;
    std::vector<catype_t> unified_catypes;

    std::vector<int3> ngbrs_14;

    std::vector<int> molecules;
    std::vector<cgrp_t> charge_groups;

    topo_t topo = {};

    /* =============================================
     * == FROM FEP FILE
     * =============================================
     */

    int n_lambdas = 0;

    int n_qangcouples = 0;
    int n_qangles = 0;
    int n_qbonds = 0;
    int n_qcangles = 0;
    int n_qcatypes = 0;
    int n_qcbonds = 0;
    int n_qcimpropers = 0;
    int n_qctorsions = 0;
    int n_qelscales = 0;
    int n_qexclpairs = 0;
    int n_qimprcouples = 0;
    int n_qimpropers = 0;
    int n_qoffdiags = 0;
    int n_qshakes = 0;
    int n_qsoftpairs = 0;
    int n_qsoftcores = 0;
    int n_qtorcouples = 0;
    int n_qtorsions = 0;

    std::vector<q_angcouple_t> q_angcouples;
    std::vector<int> q_atoms;
    std::vector<cangle_t> q_cangles;
    std::vector<catype_t> q_catypes;
    std::vector<cbond_t> q_cbonds;
    std::vector<q_cimproper_t> q_cimpropers;
    std::vector<ctorsion_t> q_ctorsions;
    std::vector<q_offdiag_t> q_offdiags;
    std::vector<q_imprcouple_t> q_imprcouples;
    std::vector<q_softpair_t> q_softpairs;
    std::vector<q_torcouple_t> q_torcouples;

    std::vector<angle_t> q_angles;
    std::vector<atype_t> q_atypes;
    std::vector<bond_t> q_bonds;
    std::vector<ccharge_t> q_charges;
    std::unique_ptr<HostDeviceBuffer<q_elscale_t>> q_elscales;
    std::vector<q_exclpair_t> q_exclpairs;
    std::vector<q_improper_t> q_impropers;
    std::vector<q_shake_t> q_shakes;
    std::vector<q_softcore_t> q_softcores;
    std::vector<torsion_t> q_torsions;

    /* =============================================
     * == RESTRAINTS
     * =============================================
     */

    /* =============================================
     * == SHELLS / SOLVENT
     * =============================================
     */

    double crgQtot = 0.0;
    double Dwmz = 0.0;
    double awmz = 0.0;

    std::vector<double> theta;
    std::vector<double> theta0;
    std::vector<double> tdum;
    int n_max_inshell = 0;
    int n_shells = 0;
    std::vector<std::vector<int>> list_sh;
    std::vector<std::vector<int>> nsort;

    /* =============================================
     * == SHAKE
     * =============================================
     */


    /* =============================================
     * == CALCULATED IN THE INTEGRATION
     * =============================================
     */

    std::vector<int> p_atoms;

    energy_t E_total = {};
    std::vector<energy_t> EQ_total;

    E_bonded_t E_bond_p = {};
    E_bonded_t E_bond_w = {};
    E_bonded_t E_bond_q = {};
    std::vector<E_bonded_t> EQ_bond;

    E_nonbonded_t E_nonbond_pp = {};
    E_nonbonded_t E_nonbond_pw = {};
    E_nonbonded_t E_nonbond_ww = {};
    E_nonbonded_t E_nonbond_qx = {};
    std::vector<E_nonbonded_t> EQ_nonbond_qq;
    std::vector<E_nonbonded_t> EQ_nonbond_qp;
    std::vector<E_nonbonded_t> EQ_nonbond_qw;
    std::vector<E_nonbonded_t> EQ_nonbond_qx;

    E_restraint_t E_restraint = {};

    double Temp = 0.0;
    double Tfree = 0.0;
    double Texcl = 0.0;
    double Tscale = 1.0;
    double A_O = 0.0;
    double A_OO = 0.0;
    double B_O = 0.0;
    double B_OO = 0.0;
    double crg_ow = 0.0;
    double crg_hw = 0.0;
    double mu_w = 0.0;

    /* =============================================
     * == ENERGY & TEMPERATURE
     * =============================================
     */

    bool separate_scaling = false;
    double Ndegf = 0.0;
    double Ndegfree = 0.0;
    double Ndegf_solvent = 0.0;
    double Ndegf_solute = 0.0;
    double Ndegfree_solvent = 0.0;
    double Ndegfree_solute = 0.0;
    double Tscale_solute = 0.0;
    double Tscale_solvent = 0.0;

    int n_parameter_states() const {
        return n_lambdas > 0 ? n_lambdas : 1;
    }

    int unified_parameter_code(int atom_idx, int state) const {
        const int qi = atom_to_qi[atom_idx];
        if (qi == -1 || state == 0) {
            return atom_idx + 1;
        }
        return n_atoms + (state - 1) * n_qatoms + qi + 1;
    }

    int unified_charge_code(int atom_idx, int state) const {
        return unified_parameter_code(atom_idx, state);
    }

    const ccharge_t& unified_ccharge_by_code(int code) const {
        return unified_ccharges[code - 1];
    }

    const ccharge_t& unified_ccharge(int atom_idx, int state) const {
        return unified_ccharge_by_code(unified_charge_code(atom_idx, state));
    }

    int unified_atype_code(int atom_idx, int state) const {
        return unified_parameter_code(atom_idx, state);
    }

    const catype_t& unified_catype_by_code(int code) const {
        return unified_catypes[code - 1];
    }

    const catype_t& unified_catype(int atom_idx, int state) const {
        return unified_catype_by_code(unified_atype_code(atom_idx, state));
    }

   private:
    Context() = default;
    ~Context() {}

    Context(const Context&) = delete;
    Context& operator=(const Context&) = delete;
};
