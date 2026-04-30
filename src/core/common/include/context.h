#pragma once

#include <vector_types.h>

#include <array>
#include <map>
#include <memory>
#include <string>
#include <vector>

#include "common/include/md_types.h"
#include "common/include/nonbonded_14_mode.h"
#include "common/include/vdw_rules.h"
#include "host_device_buffer.h"

class Context {
   public:
    Context();
    ~Context();

    static Context& instance();
    /*
    Config
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
    topo_t topo;
    int n_excluded = 0;
    int n_lambdas = 0;
    bool separate_scaling = false;
    charge_group_config_t charge_group_config;

    /*
    Atoms move variables
    */
    std::unique_ptr<HostDeviceBuffer<coord_t>> coords;
    std::unique_ptr<HostDeviceBuffer<vel_t>> velocities;
    std::unique_ptr<HostDeviceBuffer<dvel_t>> dvelocities;


    /*
    Bond forces
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

    int n_ngbrs14 = 0;
    std::unique_ptr<HostDeviceBuffer<int3>> ngbrs_14;
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

    std::vector<int> atom_to_qi;
    std::unique_ptr<HostDeviceBuffer<ccharge_t>> unified_ccharges;
    std::unique_ptr<HostDeviceBuffer<catype_t>> unified_catypes;

    std::vector<int> q_atoms;
    std::vector<int> p_atoms;
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
    std::vector<int> molecules;

    /*
    Water
    */
    std::unique_ptr<HostDeviceBuffer<shell_t>> wshells;
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


    /*
    FEP
    */
    std::unique_ptr<HostDeviceBuffer<double>> lambdas; // Actually length is only 2..

    /*
    Energy
    */

    std::unique_ptr<HostDeviceBuffer<E_restraint_t>> EQ_restraint;
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



    /*
    Pre compute Info for non bonded calculation
    */

    std::unique_ptr<HostDeviceBuffer<int>> p_atoms_list;
    std::unique_ptr<HostDeviceBuffer<int>> w_atoms_list;
    std::unique_ptr<HostDeviceBuffer<int>> q_atoms_list;
    std::unique_ptr<HostDeviceBuffer<real_t>> charge_pair_products;
    std::unique_ptr<HostDeviceBuffer<int>> p_charge_types;
    std::unique_ptr<HostDeviceBuffer<int>> w_charge_types;
    std::unique_ptr<HostDeviceBuffer<int>> q_charge_types;

    std::unique_ptr<HostDeviceBuffer<vdw_pair_param_t>> catype_pair_params;
    std::unique_ptr<HostDeviceBuffer<int>> p_catype_types;
    std::unique_ptr<HostDeviceBuffer<int>> w_catype_types;
    std::unique_ptr<HostDeviceBuffer<int>> q_catype_types;

    int n_charge_types = 0;
    int zero_charge_type = -1;
    int n_catype_types = 0;
    int zero_catype_type = -1;

    /*
    Temperature
    */

    double Temp = 0.0;
    double Tfree = 0.0;
    double Ndegf = 0.0;
    double Ndegfree = 0.0;

    double Tscale_solute = 0.0;
    double Tscale_solvent = 0.0;
    /*
    Info for FEP
    */

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
        return unified_ccharges->cpu_data_p[code - 1];
    }

    const ccharge_t& unified_ccharge(int atom_idx, int state) const {
        return unified_ccharge_by_code(unified_charge_code(atom_idx, state));
    }

    int unified_atype_code(int atom_idx, int state) const {
        return unified_parameter_code(atom_idx, state);
    }

    const catype_t& unified_catype_by_code(int code) const {
        return unified_catypes->cpu_data_p[code - 1];
    }

    const catype_t& unified_catype(int atom_idx, int state) const {
        return unified_catype_by_code(unified_atype_code(atom_idx, state));
    }

    void cuda_reset_energies();

    void init();
    void preprocess_data();


   private:
    static Context* current_;

    Context(const Context&) = delete;
    Context& operator=(const Context&) = delete;
    Context(Context&&) = delete;
    Context& operator=(Context&&) = delete;

    void init_data_from_files();
    void gpu_data_upload();

};
