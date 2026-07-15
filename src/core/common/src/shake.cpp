#include "shake.h"

#include <set>
#include <vector>

void Shake::init(Context& ctx, const ParseResult& parsed) {
    build_constraints(ctx, parsed);
    init_backend(ctx);
}

void Shake::build_constraints(Context& ctx, const ParseResult& parsed) {
    /*
     */
    data_.n_constraints = 0;
    std::vector<ShakeBond> shake_bonds;
    std::vector<int> mol_n_shakes(ctx.n_molecules());

    int current_mol = 0;

    const auto& bonds = parsed.bonds;
    const auto& cbonds = parsed.cbonds;
    const auto& heavy = ctx.heavy->cpu_data_p;
    std::set<std::pair<int, int>> constrained_pairs;

    for (int bi = 0; bi < bonds.size(); bi++) {
        int ai = bonds[bi].ai - 1;
        int aj = bonds[bi].aj - 1;

        // get the current mol index
        while (current_mol + 1 < ctx.n_molecules() && ai + 1 >= ctx.molecules[current_mol + 1]) {
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

            auto pair = std::minmax(ai, aj);
            constrained_pairs.insert(pair); 
        }
    }
    // upload
    data_.shake_bonds = HostDeviceBuffer<ShakeBond>::from_vector(shake_bonds, ctx.command_info.requested_gpu);
    data_.mol_n_shakes = HostDeviceBuffer<int>::from_vector(mol_n_shakes, ctx.command_info.requested_gpu);
}

