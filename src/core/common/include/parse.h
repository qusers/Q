#pragma once

#include <string>
#include <utility>
#include <vector>

#include "command_parser.h"
#include "md_types.h"

struct NativeOutputConfig {
    std::string final_file;
    std::string trajectory_file;
    std::string energy_file;
    std::string trajectory_atoms;
    std::string topology_file;
};

class ParseResult {
   public:
    bool fresh_start = false;
    NativeOutputConfig native_output;

    md_t md;
    topo_t topo;

    int n_atoms_solute = 0;
    int n_bonds_solute = 0;
    int n_angles_solute = 0;
    int n_torsions_solute = 0;
    int n_impropers_solute = 0;

    std::vector<double> lambdas;
    std::vector<coord_t> coords_init;
    std::vector<coord_t> coords;
    std::vector<vel_t> velocities;
    std::vector<double> restart_theta_corr;

    std::vector<bond_t> bonds;
    std::vector<cbond_t> cbonds;

    std::vector<angle_t> angles;
    std::vector<cangle_t> cangles;

    std::vector<torsion_t> torsions;
    std::vector<ctorsion_t> ctorsions;

    std::vector<improper_t> impropers;
    std::vector<cimproper_t> cimpropers;

    std::vector<restrpos_t> restrspos;
    std::vector<restrang_t> restrangs;
    std::vector<restrdis_t> restrdists;
    std::vector<restrseq_t> restrseqs;
    std::vector<restrwall_t> restrwalls;

    std::vector<charge_t> charges;
    std::vector<ccharge_t> ccharges;

    std::vector<atype_t> atypes;
    std::vector<catype_t> catypes;

    std::vector<bool> excluded;

    std::vector<std::pair<int, int>> ngbrs14;
    std::vector<std::pair<int, int>> ngbrs14_long;
    std::vector<std::pair<int, int>> ngbrs23;
    std::vector<std::pair<int, int>> ngbrs23_long;
    std::vector<int> molecules;
    charge_group_config_t charge_groups;

    std::vector<int> q_atoms;
    std::vector<q_angcouple_t> q_angcouples;

    std::vector<atype_t> q_atypes;
    std::vector<catype_t> q_catypes;

    std::vector<ccharge_t> q_charges;

    std::vector<angle_t> q_angles;
    std::vector<cangle_t> q_cangles;

    std::vector<bond_t> q_bonds;
    std::vector<cbond_t> q_cbonds;

    std::vector<q_improper_t> q_impropers;
    std::vector<q_cimproper_t> q_cimpropers;

    std::vector<torsion_t> q_torsions;
    std::vector<ctorsion_t> q_ctorsions;

    std::vector<q_offdiag_t> q_offdiags;

    std::vector<q_imprcouple_t> q_imprcouples;

    std::vector<q_softpair_t> q_softpairs;

    std::vector<q_torcouple_t> q_torcouples;

    std::vector<q_elscale_t> q_elscales;
    std::vector<q_exclpair_t> q_exclpairs;
    std::vector<q_shake_t> q_shakes;
    std::vector<q_softcore_t> q_softcores;
    bool softcore_use_max_potential = false;
};

class BaseParser {
   public:
    explicit BaseParser(const std::string& file_path) : file_path(file_path) {}
    virtual ~BaseParser() = default;

    ParseResult parse() {
        result = ParseResult{};

        parse_md();
        parse_lambdas();
        parse_topo();

        parse_coords_init();
        parse_coords();
        parse_velocities();

        parse_bonds();
        parse_cbonds();

        parse_angles();
        parse_cangles();

        parse_torsions();
        parse_ctorsions();

        parse_impropers();
        parse_cimpropers();

        parse_restrspos();
        parse_restrangs();
        parse_restrdists();
        parse_restrseqs();
        parse_restrwalls();

        parse_charges();
        parse_ccharges();

        parse_atypes();
        parse_catypes();

        parse_heavy();
        parse_excluded();

        parse_ngbrs14();
        parse_ngbrs14_long();
        parse_ngbrs23();
        parse_ngbrs23_long();

        parse_molecules();

        parse_charge_groups();

        parse_q_atoms();
        parse_q_angcouples();

        parse_q_atypes();
        parse_q_catypes();

        parse_q_charges();

        parse_q_angles();
        parse_q_cangles();

        parse_q_bonds();
        parse_q_cbonds();

        parse_q_impropers();
        parse_q_cimpropers();

        parse_q_torsions();
        parse_q_ctorsions();

        parse_q_offdiags();

        parse_q_imprcouples();

        parse_q_softpairs();

        parse_q_torcouples();

        parse_q_elscales();

        parse_q_exclpairs();

        parse_q_shakes();

        parse_q_softcores();

        return result;
    }

   protected:
    std::string file_path;
    ParseResult result;

    virtual void parse_md() = 0;
    virtual void parse_lambdas() = 0;
    virtual void parse_topo() = 0;

    virtual void parse_coords_init() = 0;
    virtual void parse_coords() = 0;
    virtual void parse_velocities() = 0;

    virtual void parse_bonds() = 0;
    virtual void parse_cbonds() = 0;

    virtual void parse_angles() = 0;
    virtual void parse_cangles() = 0;

    virtual void parse_torsions() = 0;
    virtual void parse_ctorsions() = 0;

    virtual void parse_impropers() = 0;
    virtual void parse_cimpropers() = 0;

    virtual void parse_restrspos() = 0;
    virtual void parse_restrangs() = 0;
    virtual void parse_restrdists() = 0;
    virtual void parse_restrseqs() = 0;
    virtual void parse_restrwalls() = 0;

    virtual void parse_charges() = 0;
    virtual void parse_ccharges() = 0;

    virtual void parse_atypes() = 0;
    virtual void parse_catypes() = 0;

    virtual void parse_heavy() = 0;
    virtual void parse_excluded() = 0;

    virtual void parse_ngbrs14() = 0;
    virtual void parse_ngbrs14_long() = 0;
    virtual void parse_ngbrs23() = 0;
    virtual void parse_ngbrs23_long() = 0;

    virtual void parse_molecules() = 0;

    virtual void parse_charge_groups() = 0;

    virtual void parse_q_atoms() = 0;

    virtual void parse_q_angcouples() = 0;

    virtual void parse_q_atypes() = 0;
    virtual void parse_q_catypes() = 0;

    virtual void parse_q_charges() = 0;

    virtual void parse_q_angles() = 0;
    virtual void parse_q_cangles() = 0;

    virtual void parse_q_bonds() = 0;
    virtual void parse_q_cbonds() = 0;

    virtual void parse_q_impropers() = 0;
    virtual void parse_q_cimpropers() = 0;

    virtual void parse_q_torsions() = 0;
    virtual void parse_q_ctorsions() = 0;

    virtual void parse_q_offdiags() = 0;

    virtual void parse_q_imprcouples() = 0;

    virtual void parse_q_softpairs() = 0;

    virtual void parse_q_torcouples() = 0;

    virtual void parse_q_elscales() = 0;

    virtual void parse_q_exclpairs() = 0;

    virtual void parse_q_shakes() = 0;

    virtual void parse_q_softcores() = 0;
};
