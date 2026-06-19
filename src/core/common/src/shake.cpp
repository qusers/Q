#include "shake.h"

#include <set>
#include <vector>

void Shake::init(Context& ctx) {
    build_constraints(ctx);
    exclude_constrained_bonded_terms(ctx);
    init_backend(ctx);
}

void Shake::build_constraints(Context& ctx) {
    /*
     */
    data_.n_constraints = 0;
    std::vector<ShakeBond> shake_bonds;
    std::vector<int> mol_n_shakes(ctx.n_molecules);

    int n_bonds = ctx.n_bonds;
    int current_mol = 0;

    const auto& bonds = ctx.bonds->cpu_data_p;
    const auto& cbonds = ctx.cbonds->cpu_data_p;
    const auto& heavy = ctx.heavy->cpu_data_p;

    for (int bi = 0; bi < n_bonds; bi++) {
        int ai = bonds[bi].ai - 1;
        int aj = bonds[bi].aj - 1;

        // get the current mol index
        while (current_mol + 1 < ctx.n_molecules && ai + 1 >= ctx.molecules[current_mol + 1]) {
            current_mol += 1;
        }

        if ((ctx.md.shake_hydrogens && (!heavy[ai] || !heavy[aj])) ||
            (ctx.md.shake_solute && ai + 1 <= ctx.n_atoms_solute) ||
            (ctx.md.shake_solvent && ai + 1 > ctx.n_atoms_solute)) {
            data_.n_constraints++;
            double dist = cbonds[bonds[bi].code - 1].b0;
            double dist2 = dist * dist;
            shake_bonds.emplace_back(ShakeBond{ai + 1, aj + 1, dist2});
            mol_n_shakes[current_mol]++;
        }
    }
    // upload
    data_.shake_bonds = HostDeviceBuffer<ShakeBond>::from_vector(shake_bonds, ctx.command_info.requested_gpu);
    data_.mol_n_shakes = HostDeviceBuffer<int>::from_vector(mol_n_shakes, ctx.command_info.requested_gpu);
}

void Shake::exclude_constrained_bonded_terms(Context& ctx) {
    /*
     * SHAKE constraints replace the corresponding bonded force terms.
     * The rule belongs to SHAKE, but the active bond/angle lists are owned
     * by Context, so this initialization pass intentionally rewrites them.
     */
    if (data_.n_constraints == 0) return;

    std::set<std::pair<int, int>> S;
    auto& shake_bonds = data_.shake_bonds->cpu_data_p;
    for (int i = 0; i < data_.n_constraints; i++) {
        S.insert({shake_bonds[i].ai, shake_bonds[i].aj});
    }

    int excluded = 0;
    int solute_excluded = 0;
    std::vector<bond_t> new_bonds;
    const auto& bonds = ctx.bonds->cpu_data_p;

    for (int i = 0; i < ctx.n_bonds; i++) {
        if (S.find({bonds[i].ai, bonds[i].aj}) != S.end() || S.find({bonds[i].aj, bonds[i].ai}) != S.end()) {
            excluded++;
            if (i < ctx.n_bonds_solute) solute_excluded++;
        } else {
            new_bonds.emplace_back(bonds[i]);
        }
    }
    ctx.bonds = HostDeviceBuffer<bond_t>::from_vector(new_bonds, ctx.command_info.requested_gpu);
    ctx.n_bonds -= excluded;
    ctx.n_bonds_solute -= solute_excluded;

    excluded = 0;
    solute_excluded = 0;
    std::vector<angle_t> new_angles;
    const auto& angles = ctx.angles->cpu_data_p;
    for (int i = 0; i < ctx.n_angles; i++) {
        if (S.find({angles[i].ai, angles[i].ak}) != S.end() || S.find({angles[i].ak, angles[i].ai}) != S.end()) {
            excluded++;
            if (i < ctx.n_angles_solute) solute_excluded++;
        } else {
            new_angles.emplace_back(angles[i]);
        }
    }

    ctx.angles = HostDeviceBuffer<angle_t>::from_vector(new_angles, ctx.command_info.requested_gpu);
    ctx.n_angles -= excluded;
    ctx.n_angles_solute -= solute_excluded;
}
