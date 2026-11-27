#include "temperature.h"

#include "constants.h"

void Temperature::init(Context& ctx, const ConstraintForce& constraint_force) {
    double excl_shake = 0;
    auto* excluded = ctx.excluded->cpu_data_p;
    auto* constraint_bonds = constraint_force.data().constraint_bonds->cpu_data_p;
    int n_constraints = constraint_force.data().n_constraints;

    int mol = 0, n_solute_shake_constraints = 0;
    for (int i = 0; i < n_constraints; i++) {
        const int ai = constraint_bonds[i].ai - 1;
        const int aj = constraint_bonds[i].aj - 1;
        if (excluded[ai]) excl_shake += 0.5;
        if (excluded[aj]) excl_shake += 0.5;
        while (mol + 1 < ctx.n_molecules() && ai + 1 >= ctx.molecules[mol + 1]) mol += 1;
        if (mol < ctx.n_molecules() - ctx.n_waters()) n_solute_shake_constraints++;
    }

    data_.Ndegf = 3 * ctx.n_atoms - n_constraints;
    data_.Ndegfree = data_.Ndegf - 3 * ctx.n_excluded + excl_shake;
    data_.Ndegf_solvent = 3 * (ctx.n_atoms - ctx.n_atoms_solute) - (n_constraints - n_solute_shake_constraints);
    data_.Ndegf_solute = data_.Ndegf - data_.Ndegf_solvent;
    data_.Ndegfree_solvent = 3 * (ctx.n_atoms - ctx.n_atoms_solute) - (n_constraints - n_solute_shake_constraints);  // No frozen atoms in the solvent
    data_.Ndegfree_solute = data_.Ndegfree - data_.Ndegfree_solvent;

    if (data_.Ndegfree_solvent * data_.Ndegfree_solute == 0) {
        data_.separate_scaling = false;
    } else {
        data_.separate_scaling = ctx.md.separate_scaling;
    }
    data_.tau_T = time_unit * ctx.md.bath_coupling;
    data_.results = std::make_unique<HostDeviceBuffer<double>>(
        N_TRESULT, /*cpu*/ true, /*gpu*/ ctx.command_info.requested_gpu);

    ctx.Ndegfree = data_.Ndegfree;
    init_backend(ctx);
}

void Temperature::sync_for_output(Context& ctx) {
    const double* r = data_.results->cpu_data_p;
    ctx.Temp = r[R_TEMP];
    ctx.energy.data().Ukin = r[R_UKIN];
}