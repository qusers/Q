#pragma once

#include <memory>
#include <string>
#include <vector>

#include "parse.h"

class InpParser : public BaseParser {
   public:
    explicit InpParser(const std::string& inp_file);
    ~InpParser() override;

   protected:
    void parse_md() override;
    void parse_lambdas() override;
    void parse_topo() override;

    void parse_coords_init() override;
    void parse_coords() override;
    void parse_velocities() override;

    void parse_bonds() override;
    void parse_cbonds() override;

    void parse_angles() override;
    void parse_cangles() override;

    void parse_torsions() override;
    void parse_ctorsions() override;

    void parse_impropers() override;
    void parse_cimpropers() override;

    void parse_restrspos() override;
    void parse_restrangs() override;
    void parse_restrdists() override;
    void parse_restrseqs() override;
    void parse_restrwalls() override;

    void parse_charges() override;
    void parse_ccharges() override;

    void parse_atypes() override;
    void parse_catypes() override;

    void parse_heavy() override;
    void parse_excluded() override;

    void parse_ngbrs14() override;
    void parse_ngbrs14_long() override;
    void parse_ngbrs23() override;
    void parse_ngbrs23_long() override;

    void parse_molecules() override;

    void parse_charge_groups() override;

    void parse_q_atoms() override;

    void parse_q_angcouples() override;

    void parse_q_atypes() override;
    void parse_q_catypes() override;

    void parse_q_charges() override;

    void parse_q_angles() override;
    void parse_q_cangles() override;

    void parse_q_bonds() override;
    void parse_q_cbonds() override;

    void parse_q_impropers() override;
    void parse_q_cimpropers() override;

    void parse_q_torsions() override;
    void parse_q_ctorsions() override;

    void parse_q_offdiags() override;

    void parse_q_imprcouples() override;

    void parse_q_softpairs() override;

    void parse_q_torcouples() override;

    void parse_q_elscales() override;

    void parse_q_exclpairs() override;

    void parse_q_shakes() override;

    void parse_q_softcores() override;

   private:
    struct MDInput;
    struct TopData;
    struct FepData;

    void ensure_md_input();
    void ensure_topology();
    void ensure_fep();
    void ensure_run_start_vectors();

    std::unique_ptr<MDInput> input_;
    std::unique_ptr<TopData> top_;
    std::unique_ptr<FepData> fep_;
    std::vector<std::vector<std::string>> run_coords_;
    std::vector<std::vector<std::string>> run_velocities_;

    std::string input_dir_;
    std::string topology_file_;
    std::string fep_file_;
    std::string restart_file_;
    bool run_start_ready_ = false;
};
