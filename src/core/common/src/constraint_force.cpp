#include "constraint_force.h"

#include <set>
#include <vector>

void ConstraintForce::init(Context& ctx, const ParseResult& parsed) {
    build_constraints(ctx, parsed);
    init_backend(ctx);
}

void ConstraintForce::build_constraints(Context& ctx, const ParseResult& parsed) {
    /*
     */
    data_.n_constraints = 0;
    std::vector<ConstraintBond> constraint_bonds;
    std::vector<int> mol_n_constraints(ctx.n_molecules());

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

        if ((ctx.md.constraint_hydrogens && (!heavy[ai] || !heavy[aj])) ||
            (ctx.md.constraint_solute && ai + 1 <= ctx.n_atoms_solute) ||
            (ctx.md.constraint_solvent && ai + 1 > ctx.n_atoms_solute)) {
            data_.n_constraints++;
            double dist = cbonds[bonds[bi].code - 1].b0;
            double dist2 = dist * dist;
            constraint_bonds.emplace_back(ConstraintBond{ai + 1, aj + 1, dist2});
            mol_n_constraints[current_mol]++;

            auto pair = std::minmax(bonds[bi].ai, bonds[bi].aj);
            constrained_pairs.insert(pair);
        }
    }
    // upload
    data_.constraint_bonds = HostDeviceBuffer<ConstraintBond>::from_vector(constraint_bonds, ctx.command_info.requested_gpu);
    data_.mol_n_constraints = HostDeviceBuffer<int>::from_vector(mol_n_constraints, ctx.command_info.requested_gpu);
    data_.constrained_pairs = std::move(constrained_pairs);
}
