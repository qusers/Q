#pragma once

#include <vector_types.h>

#include <array>
#include <map>
#include <memory>
#include <string>
#include <vector>

#include "command_parser.h"
#include "common/include/md_types.h"
#include "common/include/parse.h"
#include "common/include/vdw_rules.h"
#include "host_device_buffer.h"
#include "energy.h"

class Shake;

class Context {
   public:
    Context();
    ~Context();

    static Context& instance();
    /*
    Config
    */
    CommandInfo command_info;

    bool fresh_start;
    int n_atoms;         // the total number of atoms
    int n_atoms_solute;  // the total number of solute number, in our system [0, n_atoms_solute) are solute, [n_atoms_solute, n_atoms) are water atoms
    int n_patoms;
    int n_qatoms;
    int n_waters = 0;
    int n_molecules;
    double dt = 0.0;
    double tau_T = 0.0;
    md_t md;
    topo_t topo;
    NativeOutputConfig native_output;
    int n_excluded;
    int n_lambdas;
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
    int n_angles;
    int n_angles_solute;
    int n_cangles;
    std::vector<angle_t> angles;
    std::vector<cangle_t> cangles;

    int n_bonds;
    int n_bonds_solute;
    int n_cbonds;
    std::vector<bond_t> bonds;
    std::vector<cbond_t> cbonds;
    int n_torsions;
    int n_torsions_solute;
    int n_ctorsions;
    std::vector<torsion_t> torsions;
    std::vector<ctorsion_t> ctorsions;
    int n_impropers;
    int n_impropers_solute;
    int n_cimpropers;
    std::vector<improper_t> impropers;
    std::vector<cimproper_t> cimpropers;

    int n_restrspos;
    std::vector<restrpos_t> restrspos;

    int n_restrangs;
    std::vector<restrang_t> restrangs;

    int n_restrdists;
    std::vector<restrdis_t> restrdists;

    int n_restrseqs;
    std::vector<restrseq_t> restrseqs;
    int n_restrwalls;
    std::vector<restrwall_t> restrwalls;

    /*
    Atom Info
    */
    int n_charges;
    int n_ccharges;
    std::unique_ptr<HostDeviceBuffer<charge_t>> charges;
    std::unique_ptr<HostDeviceBuffer<ccharge_t>> ccharges;
    int n_atypes;
    int n_catypes;
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

    std::unique_ptr<HostDeviceBuffer<coord_t>> xcoords;  // todo: It's just a temporary variables...
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
    std::vector<double> restart_theta_corr;
    int n_max_inshell;
    int n_shells;
    std::vector<std::vector<int>> list_sh;
    std::vector<std::vector<int>> nsort;

    /*
    FEP
    */
    std::unique_ptr<HostDeviceBuffer<double>> lambdas;  // Actually length is only 2..

    /*
    Energy
    */
    EnergyBuffer energy;


    /*
    Temperature
    */

    double Temp = 0.0;
    double Tfree = 0.0;
    int Ndegf = 0;
    int Ndegfree = 0;
    int Ndegf_solute = 0;
    int Ndegfree_solute = 0;
    int Ndegf_solvent = 0;
    int Ndegfree_solvent = 0;

    double Tscale_solute = 0.0;
    double Tscale_solvent = 0.0;
    /*
    Info for FEP
    */

    int n_qangcouples;
    int n_qangles;
    int n_qbonds;
    int n_qcangles;
    int n_qcatypes;
    int n_qcbonds;
    int n_qcimpropers;
    int n_qctorsions;
    int n_qelscales;
    int n_qexclpairs;
    int n_qimprcouples;
    int n_qimpropers;
    int n_qoffdiags;
    int n_qshakes;
    int n_qsoftpairs;
    int n_qsoftcores;
    bool softcore_use_max_potential = false;
    int n_qtorcouples;
    int n_qtorsions;

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
    void preprocess_data(Shake &shake);

   private:
    static Context* current_;

    Context(const Context&) = delete;
    Context& operator=(const Context&) = delete;
    Context(Context&&) = delete;
    Context& operator=(Context&&) = delete;

    void init_data_from_files();
};
