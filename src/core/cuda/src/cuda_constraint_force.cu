#include "cuda_constraint_force.cuh"
#include "cuda_lincs.cuh"
#include "cuda_settle.cuh"
#include "cuda_shake.cuh"
#include "helpers.h"

namespace {
std::unique_ptr<ConstraintForce> get_constraint_force_by_name(const std::string& name) {
    if (name == "shake") {
        return std::make_unique<CudaShake>();
    } else if (name == "lincs") {
        return std::make_unique<CudaLincs>();
    } else if (name == "settle") {
        return std::make_unique<CudaSettle>();
    }
    return nullptr;
}

}  // namespace

void CudaConstraintForce::init_backend(Context& ctx) {
    const std::string& solute_constraint_algorithm = ctx.md.solute_constraint_algorithm;
    const std::string& solvent_constraint_algorithm = ctx.md.solvent_constraint_algorithm;

    /*
    two cases: if they are same? we should use on kernel to calculate them.
    If they are not same, use two different kernels.
    Position	Allowed methods
    Solute	    shake, lincs
    Solvent	    shake, lincs, settle
    */

    if (solute_constraint_algorithm != "shake" && solute_constraint_algorithm != "lincs") {
        fatal(
            "[Constraint Force] Invalid solute constraint algorithm '" +
            solute_constraint_algorithm + "'; expected shake or lincs.");
    }

    if (solvent_constraint_algorithm != "shake" &&
        solvent_constraint_algorithm != "lincs" &&
        solvent_constraint_algorithm != "settle") {
        fatal(
            "[Constraint Force] Invalid solvent constraint algorithm '" +
            solvent_constraint_algorithm + "'; expected shake, lincs, or settle.");
    }

    printf("Solute constraint algorithm is %s. Solvent constraint algorithm is %s.\n",
           solute_constraint_algorithm.c_str(), solvent_constraint_algorithm.c_str());

    const int bond_num = data_.n_constraints;
    if (bond_num == 0) return;

    std::vector<ConstraintBond> constraint_bonds(data_.constraint_bonds->cpu_data_p, data_.constraint_bonds->cpu_data_p + bond_num);

    if (solute_constraint_algorithm == solvent_constraint_algorithm) {
        // shake or lincs
        common_constraint_force_ = get_constraint_force_by_name(solute_constraint_algorithm);
        common_constraint_force_->init_from_bonds(ctx, constraint_bonds);

    } else {
        std::vector<ConstraintBond> solute_bonds, solvent_bonds;

        for (int i = 0; i < bond_num; i++) {
            const int ai = constraint_bonds[i].ai - 1;
            const int aj = constraint_bonds[i].aj - 1;
            const bool ai_is_solute = ai < ctx.n_atoms_solute;
            const bool aj_is_solute = aj < ctx.n_atoms_solute;

            if (ai_is_solute != aj_is_solute) {
                fatal(
                    "[Constraint Force] Constraint between atoms " +
                    std::to_string(ai + 1) + " and " + std::to_string(aj + 1) +
                    " crosses the solute-solvent boundary.");
            }

            if (!ai_is_solute) {
                // water
                solvent_bonds.emplace_back(constraint_bonds[i]);
            } else {
                solute_bonds.emplace_back(constraint_bonds[i]);
            }
        }

        if (!solute_bonds.empty()) {
            solute_constraint_force_ = get_constraint_force_by_name(solute_constraint_algorithm);
            solute_constraint_force_->init_from_bonds(ctx, solute_bonds);
        }
        if (!solvent_bonds.empty()) {
            solvent_constraint_force_ = get_constraint_force_by_name(solvent_constraint_algorithm);
            solvent_constraint_force_->init_from_bonds(ctx, solvent_bonds);
        }
    }
}

void CudaConstraintForce::apply(Context& ctx, HostDeviceBuffer<coord_t>& xcoords) {
    if (common_constraint_force_) {
        common_constraint_force_->apply(ctx, xcoords);
    } else {
        if (solute_constraint_force_) solute_constraint_force_->apply(ctx, xcoords);
        if (solvent_constraint_force_) solvent_constraint_force_->apply(ctx, xcoords);
    }
}

void CudaConstraintForce::initial_constraint(Context& ctx) {
    if (common_constraint_force_) {
        common_constraint_force_->initial_constraint(ctx);
    } else {
        if (solute_constraint_force_) solute_constraint_force_->initial_constraint(ctx);
        if (solvent_constraint_force_) solvent_constraint_force_->initial_constraint(ctx);
    }
}

void CudaConstraintForce::cleanup() {}

void CudaConstraintForce::init_from_bonds(Context& ctx, const std::vector<ConstraintBond>& bonds) {};
