#include "cuda_lincs.cuh"
#include "vector"

void CudaLincs::apply(Context& ctx, HostDeviceBuffer<coord_t>& xcoords) {
    coord_t* candidate = ctx.coords->gpu_data_p;
    const coord_t* reference = xcoords.gpu_data_p;
    apply_to(ctx, candidate, reference);
}

void CudaLincs::init_backend(Context&) {
    int bond_num = data_.n_constraints;
    std::vector<CudaLincsBond> bonds(bond_num);
    auto* constraint_bonds = data_.constraint_bonds->cpu_data_p;
    for (int i = 0; i < bond_num; i++) {
        bonds[i] = CudaLincsBond{data_.}


    }
    




}

void CudaLincs::initial_constraint(Context&) {
}

void CudaLincs::apply_to(Context& ctx, coord_t* coords, const coord_t* xcoords) {
}

void CudaLincs::cleanup() {
}