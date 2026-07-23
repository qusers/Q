#pragma once
#include <array>
#include <cstddef>
#include <map>
#include <memory>
#include <string>
#include <vector>

#include "command_parser.h"
#include "common/include/md_types.h"
#include "common/include/parse.h"
#include "common/include/vdw_rules.h"
#include "energy.h"
#include "host_device_buffer.h"

class Shake;

class Context {
   public:
    Context() = default;
    ~Context() = default;

    CommandInfo command_info;
    bool fresh_start = false;
    int n_atoms = 0;         // the total number of atoms
    int n_atoms_solute = 0;  // the total number of solute number, in our system [0, n_atoms_solute) are solute, [n_atoms_solute, n_atoms) are water atoms
    double dt = 0.0;
    int n_excluded = 0;
    md_t md;
    topo_t topo;
    std::vector<NativeOutputConfig> native_outputs;
    charge_group_config_t charge_group_config;
    std::unique_ptr<HostDeviceBuffer<coord_t>> coords;
    std::unique_ptr<HostDeviceBuffer<vel_t>> velocities;
    std::unique_ptr<HostDeviceBuffer<dvel_t>> dvelocities;
    std::unique_ptr<HostDeviceBuffer<charge_t>> charges;
    std::unique_ptr<HostDeviceBuffer<ccharge_t>> ccharges;
    std::unique_ptr<HostDeviceBuffer<atype_t>> atypes;
    std::unique_ptr<HostDeviceBuffer<catype_t>> catypes;
    std::unique_ptr<HostDeviceBuffer<bool>> heavy;
    std::unique_ptr<HostDeviceBuffer<coord_t>> coords_init;
    std::unique_ptr<HostDeviceBuffer<bool>> excluded;
    std::unique_ptr<HostDeviceBuffer<double>> winv;
    std::unique_ptr<HostDeviceBuffer<bool>> shell;
    std::vector<int> q_atoms;
    std::vector<int> p_atoms;
    std::unique_ptr<HostDeviceBuffer<int>> LJ_matrix;
    std::vector<int> molecules;
    std::vector<std::vector<shell_t>> replica_wshells;
    std::unique_ptr<HostDeviceBuffer<double>> lambdas;  // Actually length is only 2..
    EnergyBuffer energy;
    std::vector<double> replica_temperatures;

    int Ndegfree = 0;
    std::vector<catype_t> q_catypes;
    std::vector<atype_t> q_atypes;
    std::vector<ccharge_t> q_charges;
    std::vector<ParseResult> get_parse_results();
    void init(const std::vector<ParseResult>& replica_results);
    int n_lambdas() const { return lambdas->length; }
    int n_waters() const { return (n_atoms - n_atoms_solute) / 3; }
    int n_qatoms() const { return q_atoms.size(); }
    int n_molecules() const { return molecules.size(); }
    int n_patoms() const { return p_atoms.size(); }
    int n_replicates() const { return command_info.n_replicates(); }
    double temperature(int replica = 0) const {
        return replica_temperatures[replica];
    }
    size_t n_total_atoms() const { return static_cast<size_t>(n_atoms) * static_cast<size_t>(n_replicates()); }

   private:
    void set_lj_pairs(std::vector<int>& matrix, const std::vector<std::pair<int, int>>& pair, int value);

    void init_fresh_start(const ParseResult& parsed);
    void init_md(const ParseResult& parsed);
    void init_topo(const ParseResult& parsed);
    void init_native_outputs(const std::vector<ParseResult>& parsed);
    void init_const(const ParseResult& parsed);
    void init_coords_init(const ParseResult& parsed);
    void init_coords(const std::vector<ParseResult>& replica_results);
    void init_lambdas(const ParseResult& parsed);
    void init_charges(const ParseResult& parsed);
    void init_ccharges(const ParseResult& parsed);
    void init_atypes(const ParseResult& parsed);
    void init_catypes(const ParseResult& parsed);
    void init_excluded(const ParseResult& parsed);
    void init_molecules(const ParseResult& parsed);
    void init_charge_group_config(const ParseResult& parsed);
    void init_q_atoms(const ParseResult& parsed);
    void init_q_catypes(const ParseResult& parsed);
    void init_q_atypes(const ParseResult& parsed);
    void init_q_charges(const ParseResult& parsed);
    void init_lj_matrix(const ParseResult& parsed);
    void init_velocities(const std::vector<ParseResult>& replica_results);
    void init_inv_mass();
    void init_heavy();
    void init_shell();
    void init_patoms();
    void init_dvelocities();
    void init_energy();
};
