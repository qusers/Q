#pragma once
#include <fstream>
#include <iomanip>
#include <iostream>

#include "context.h"

inline void output_coord(std::ofstream& file, const coord_t& x) {
    file << x.x << ' ' << x.y << ' ' << x.z << '\n';
}

inline void output_charge_group_config(std::ofstream& file, const charge_group_config_t& x) {
    file << x.n_cgrps_solute << ' ' << x.n_cgrps_solvent << x.iuse_switch_atom << '\n';
    file << x.charge_groups.size() << '\n';
    for (int i = 0; i < x.charge_groups.size(); i++) {
        file << x.charge_groups[i].iswitch << ' ' << x.charge_groups[i].atoms.size() << '\n';
        for (int j = 0; j < x.charge_groups[i].atoms.size(); j++) {
            file << x.charge_groups[i].atoms[j] << ' ';
        }
        file << '\n';
    }
}

inline void output_coords(std::ofstream& file, const HostDeviceBuffer<coord_t>& coords) {
    file << coords.length << '\n';
    for (int i = 0; i < coords.length; i++) {
        output_coord(file, coords.cpu_data_p[i]);
    }
}

inline void output_velocity(std::ofstream& file, const vel_t& x) {
    file << x.x << ' ' << x.y << ' ' << x.z << '\n';
}

inline void output_velocities(std::ofstream& file, const HostDeviceBuffer<vel_t>& velocities) {
    file << velocities.length << '\n';
    for (int i = 0; i < velocities.length; i++) {
        output_velocity(file, velocities.cpu_data_p[i]);
    }
}

inline void output_dvelocity(std::ofstream& file, const dvel_t& x) {
    file << x.x << ' ' << x.y << ' ' << x.z << '\n';
}

inline void output_dvelocities(std::ofstream& file, const HostDeviceBuffer<dvel_t>& dvelocities) {
    file << dvelocities.length << '\n';
    for (int i = 0; i < dvelocities.length; i++) {
        output_dvelocity(file, dvelocities.cpu_data_p[i]);
    }
}

inline void output_ctx_in_file(const Context& ctx) {
    std::ofstream file("./cuda_context_debug.txt");
    file << std::fixed << std::setprecision(8);
    file << ctx.fresh_start << '\n';
    file << ctx.n_atoms << '\n';
    file << ctx.n_atoms_solute << '\n';
    file << ctx.n_patoms << '\n';
    file << ctx.n_qatoms << '\n';
    file << ctx.n_waters << '\n';
    file << ctx.n_molecules << '\n';

    // file << ctx.dt << '\n';
    // file << ctx.tau_T << '\n';
    // file << ctx.md.steps << '\n';
    // file << ctx.md.stepsize << '\n';
    // file << ctx.md.temperature << '\n';
    // file << ctx.md.bath_coupling << '\n';
    // file << ctx.md.random_seed << '\n';  // different when it's not beginning state. fortran is 0, and we are 1.
    // file << ctx.md.initial_temperature << '\n';
    // file << ctx.md.shake_solvent << '\n';
    // file << ctx.md.shake_solute << '\n';
    // file << ctx.md.shake_hydrogens << '\n';
    // file << ctx.md.lrf << '\n';
    // file << ctx.md.separate_scaling << '\n';
    // file << ctx.md.charge_groups << '\n';
    // file << ctx.md.solute_solute << '\n';
    // file << ctx.md.solvent_solvent << '\n';
    // file << ctx.md.solute_solvent << '\n';
    // file << ctx.md.q_atom << '\n';
    // file << ctx.md.shell_radius << '\n';
    // file << ctx.md.shell_force << '\n';
    // file << ctx.md.radial_force << '\n';
    // file << ctx.md.polarisation << '\n';
    // file << ctx.md.polarisation_force << '\n';
    // file << ctx.md.non_bond << '\n';
    // file << ctx.md.output << '\n';
    // file << ctx.md.energy << '\n';
    // file << ctx.md.trajectory << '\n';
    // file << ctx.topo.solvent_type << '\n';
    // file << ctx.topo.exclusion_radius << '\n';
    // file << ctx.topo.solvent_radius << '\n';
    // output_coord(file, ctx.topo.solute_center);
    // output_coord(file, ctx.topo.solvent_center);
    // file << ctx.topo.el14_scale << '\n';
    // file << ctx.topo.coulomb_constant << '\n';

    // file << ctx.n_excluded << '\n';
    // file << ctx.n_lambdas << '\n';
    // file << ctx.separate_scaling << '\n';
    // output_charge_group_config(file, ctx.charge_group_config);
    // output_coords(file, *ctx.coords);
    output_velocities(file, *ctx.velocities);
    output_dvelocities(file, *ctx.dvelocities);

    file.close();
}
inline void debug(const std::string& s) {
    static bool flag = false;
    if (!flag) {
        std::ofstream("./debug.txt", std::ios::trunc).close();
        flag = true;
    }
    std::ofstream file("./debug.txt", std::ios::app);
    file << std::fixed << std::setprecision(8);
    file << s << '\n';
}