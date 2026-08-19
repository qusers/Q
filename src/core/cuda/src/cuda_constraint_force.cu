#include "cuda_constraint_force.cuh"
#include "cuda_lincs.cuh"
#include "cuda_settle.cuh"
#include "cuda_shake.cuh"

namespace {
std::unique_ptr<ConstraintForce> get_constraint_force_by_name(std::string name) {
    if (name == "shake") {
        return std::unique_ptr<ConstraintForce>(std::make_unique<CudaShake>());
    } else if (name == "lincs") {
        return std::unique_ptr<ConstraintForce>(std::make_unique<CudaLincs>());
    } else if (name == "settle") {
        return std::unique_ptr<ConstraintForce>(std::make_unique<CudaSettle>());
    }
    return nullptr;
}

}  // namespace

void CudaConstraintForce::init_backend(Context& ctx) {
    auto solute_constraint_algorithm = ctx.md.solute_constraint_algorithm;
    auto solvent_constraint_algorithm = ctx.md.solvent_constraint_algorithm;

    // todo: read from md
    solute_constraint_algorithm = "lincs";
    solvent_constraint_algorithm = "settle";

    /*
    two cases: if they are same? we should use on kernel to calculate them.
    If they are not same, use two different kernels.
    Position	Allowed methods
    Solute	    shake, lincs
    Solvent	    shake, lincs, settle
    */

    std::set<std::string> valid_algorithm = {"shake", "lincs", "settle"};

    if (valid_algorithm.count(solute_constraint_algorithm) == 0) {
        printf("[Constraint Force] Invalid solute constraint algorithm : %s, only support shake, lincs, settle", solute_constraint_algorithm.c_str());
        std::exit(EXIT_FAILURE);
    }

    if (valid_algorithm.count(solvent_constraint_algorithm) == 0) {
        printf("[Constraint Force] Invalid solvent constraint algorithm : %s, only support shake, lincs, settle", solvent_constraint_algorithm.c_str());
        std::exit(EXIT_FAILURE);
    }

    if (solute_constraint_algorithm == "settle") {
        printf("[Constraint Force] Don't support settle for solute!");
        std::exit(EXIT_FAILURE);
    }

    int bond_num = data_.n_constraints;
    std::vector<ConstraintBond> constraint_bonds(data_.constraint_bonds->cpu_data_p, data_.constraint_bonds->cpu_data_p + bond_num);

    if (solute_constraint_algorithm == solvent_constraint_algorithm) {
        // shake or lincs
        common_constraint_force_ = get_constraint_force_by_name(solute_constraint_algorithm);
        common_constraint_force_->init_from_bonds(ctx, constraint_bonds);

    } else {
        solute_constraint_force_ = get_constraint_force_by_name(solute_constraint_algorithm);
        solvent_constraint_force_ = get_constraint_force_by_name(solvent_constraint_algorithm);

        std::vector<ConstraintBond> solute_bonds, solvent_bonds;

        for (int i = 0; i < bond_num; i++) {
            int ai = constraint_bonds[i].ai - 1, aj = constraint_bonds[i].aj - 1;
            if (ai >= ctx.n_atoms_solute) {
                // water
                solvent_bonds.emplace_back(constraint_bonds[i]);
            } else {
                solute_bonds.emplace_back(constraint_bonds[i]);
            }
        }
        solute_constraint_force_->init_from_bonds(ctx, solute_bonds);
        solvent_constraint_force_->init_from_bonds(ctx, solvent_bonds);
    }
}

void CudaConstraintForce::apply(Context& ctx, HostDeviceBuffer<coord_t>& xcoords) {
    if (common_constraint_force_) {
        common_constraint_force_->apply(ctx, xcoords);
    } else {
        solute_constraint_force_->apply(ctx, xcoords);
        solvent_constraint_force_->apply(ctx, xcoords);
    }
}

void CudaConstraintForce::initial_constraint(Context& ctx) {
    // Same idea as the CudaShake
    HostDeviceBuffer<coord_t> xcoords(ctx.n_atoms, true, true);
    auto* d_coords = ctx.coords->gpu_data_p;
    auto* d_xcoords = xcoords.gpu_data_p;

    check_cuda(cudaMemcpy(d_xcoords, d_coords, ctx.n_atoms * sizeof(coord_t), cudaMemcpyDeviceToDevice));

    apply(ctx, xcoords);

    ctx.coords->download();

    auto* coords = ctx.coords->cpu_data_p;
    auto* velocities = ctx.velocities->cpu_data_p;
    auto* host_xcoords = xcoords.cpu_data_p;
    for (int i = 0; i < ctx.n_atoms; i++) {
        const double dt = ctx.dt;
        host_xcoords[i].x = coords[i].x - dt * velocities[i].x;
        host_xcoords[i].y = coords[i].y - dt * velocities[i].y;
        host_xcoords[i].z = coords[i].z - dt * velocities[i].z;
    }

    xcoords.upload();

    apply(ctx, xcoords);

    xcoords.download();

    for (int i = 0; i < ctx.n_atoms; i++) {
        const double dt = ctx.dt;
        velocities[i].x = (coords[i].x - host_xcoords[i].x) / dt;
        velocities[i].y = (coords[i].y - host_xcoords[i].y) / dt;
        velocities[i].z = (coords[i].z - host_xcoords[i].z) / dt;
    }
    ctx.velocities->upload();
}

void CudaConstraintForce::cleanup() {}

void CudaConstraintForce::init_from_bonds(Context& ctx, const std::vector<ConstraintBond>& bonds) {};