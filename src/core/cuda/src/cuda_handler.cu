#include <iostream>
#include <vector>

#include "common/include/context.h"
#include "common/include/constants.h"
#include "common/include/precision.h"
#include "cpu/include/cpu_nonbonded_pp_force.h"
#include "cpu/include/cpu_nonbonded_pw_force.h"
#include "cpu/include/cpu_nonbonded_qp_force.h"
#include "cpu/include/cpu_nonbonded_qq_force.h"
#include "cpu/include/cpu_nonbonded_qw_force.h"
#include "cpu/include/cpu_nonbonded_ww_force.h"
#include "cpu/include/cpu_polx_water_force.h"
#include "cpu/include/cpu_temperature.h"
#include "cpu/include/cpu_q_angle_force.h"
#include "cpu/include/cpu_q_bond_force.h"
#include "cpu/include/cpu_q_improper_force.h"
#include "cpu/include/cpu_q_torsion_force.h"
#include "cpu/include/cpu_radix_water_force.h"
#include "cuda/include/cuda_angle_force.cuh"
#include "cuda/include/cuda_bond_force.cuh"
#include "cuda/include/cuda_handler.cuh"
#include "cuda/include/cuda_improper2_force.cuh"
#include "cuda/include/cuda_leapfrog.cuh"
#include "cuda/include/cuda_nonbonded_force.cuh"
#include "cuda/include/cuda_nonbonded_14_force.cuh"
#include "cuda/include/cuda_nonbonded_pp_force.cuh"
#include "cuda/include/cuda_nonbonded_pw_force.cuh"
#include "cuda/include/cuda_nonbonded_qp_force.cuh"
#include "cuda/include/cuda_nonbonded_qq_force.cuh"
#include "cuda/include/cuda_nonbonded_qw_force.cuh"
#include "cuda/include/cuda_nonbonded_ww_force.cuh"
#include "cuda/include/cuda_polx_water_force.cuh"
#include "cuda/include/cuda_pshell_force.cuh"
#include "cuda/include/cuda_radix_water_force.cuh"
#include "cuda/include/cuda_restrang_force.cuh"
#include "cuda/include/cuda_restrdis_force.cuh"
#include "cuda/include/cuda_restrpos_force.cuh"
#include "cuda/include/cuda_restrseq_force.cuh"
#include "cuda/include/cuda_restrwall_force.cuh"
#include "cuda/include/cuda_temperature.cuh"
#include "cuda/include/cuda_torsion_force.cuh"
#include "cuda/include/cuda_force_accum.cuh"
#include "cuda_utility.cuh"

namespace {

void add_cpu_forces_to_cuda_and_clear(Context& host) {
    std::vector<cuda_dvel_t> gpu_forces(host.n_atoms);
    check_cuda(cudaMemcpy(
        gpu_forces.data(),
        cuda_force_accum_buffer(host),
        sizeof(cuda_dvel_t) * gpu_forces.size(),
        cudaMemcpyDeviceToHost));

    auto* cpu_forces = host.dvelocities->cpu_data_p;
    auto add_force_component = [](force_fixed_storage_t stored, force_accum_t delta) {
        const force_fixed_t sum = force_fixed_from_storage(stored) + force_to_fixed(static_cast<double>(delta));
        return force_fixed_to_storage(sum);
    };
    for (int i = 0; i < host.n_atoms; i++) {
        auto& cpu_force = cpu_forces[i];
        gpu_forces[i].x = add_force_component(gpu_forces[i].x, cpu_force.x);
        gpu_forces[i].y = add_force_component(gpu_forces[i].y, cpu_force.y);
        gpu_forces[i].z = add_force_component(gpu_forces[i].z, cpu_force.z);
        cpu_force.x = 0;
        cpu_force.y = 0;
        cpu_force.z = 0;
    }

    check_cuda(cudaMemcpy(
        cuda_force_accum_buffer(host),
        gpu_forces.data(),
        sizeof(cuda_dvel_t) * gpu_forces.size(),
        cudaMemcpyHostToDevice));
}

void calc_q_bonded_forces_host_fallback(Context& host) {
    if (host.n_qatoms == 0 || host.n_lambdas == 0) return;
    if (host.n_qangles == 0 && host.n_qbonds == 0 && host.n_qtorsions == 0 && host.n_qimpropers == 0) return;

    for (int state = 0; state < host.n_lambdas; state++) {
        calc_qbond_forces(state);
        calc_qangle_forces(state);
        calc_qtorsion_forces(state);
        calc_qimproper_forces(state);
    }
    add_cpu_forces_to_cuda_and_clear(host);
}

void calc_q_q_nonbonded_forces_host_fallback(Context& host) {
    if (host.n_qatoms == 0 || host.n_lambdas == 0) return;

    calc_nonbonded_qq_forces();
    add_cpu_forces_to_cuda_and_clear(host);
}

}  // namespace

void CudaHandler::initialize_backend() {
    if (!initialized_) {
        init_angle_force_kernel_data();
        init_bond_force_kernel_data();
        init_improper2_force_kernel_data();
        init_leapfrog_kernel_data();
        init_nonbonded_force_kernel_data();
        init_nonbonded_14_force_kernel_data();
        init_nonbonded_pp_force_kernel_data();
        init_nonbonded_pw_force_kernel_data();
        init_nonbonded_qp_force_kernel_data();
        init_nonbonded_qq_force_kernel_data();
        init_nonbonded_qw_force_kernel_data();
        init_nonbonded_ww_force_kernel_data();
        init_polx_water_force_kernel_data();
        init_pshell_force_kernel_data();
        init_radix_water_force_kernel_data();
        init_restrang_force_kernel_data();
        init_restrdis_force_kernel_data();
        init_restrpos_force_kernel_data();
        init_restrseq_force_kernel_data();
        init_restrwall_force_kernel_data();
        init_temperature_kernel_data();
        init_torsion_force_kernel_data();
        initialized_ = true;
    }
}

void CudaHandler::shutdown() {
    if (initialized_) {
        cleanup_angle_force();
        cleanup_bond_force();
        cleanup_improper2_force();
        cleanup_leapfrog();
        cleanup_nonbonded_pp_force();
        cleanup_nonbonded_pw_force();
        cleanup_nonbonded_qp_force();
        cleanup_nonbonded_qq_force();
        cleanup_nonbonded_qw_force();
        cleanup_nonbonded_ww_force();
        cleanup_nonbonded_14_force();
        cleanup_nonbonded_force();
        cleanup_polx_water_force();
        cleanup_pshell_force();
        cleanup_radix_water_force();
        cleanup_restrang_force();
        cleanup_restrdis_force();
        cleanup_restrpos_force();
        cleanup_restrseq_force();
        cleanup_restrwall_force();
        cleanup_temperature();
        cleanup_torsion_force();
        initialized_ = false;
    }
}

void CudaHandler::calc_internal_forces(int iteration) {
    auto& host = Context::instance();
    calc_nonbonded_pp_forces_host_v2();
    if (host.topo.vdw_rule == VDW_ARITHMETIC) {
        calc_nonbonded_qp_forces_host_v2();
        calc_nonbonded_pw_forces_host_v2();
        calc_nonbonded_qw_forces_host_v2();
        calc_nonbonded_ww_forces_host_v2();
    } else {
        calc_nonbonded_pw_forces_host_v2();
        calc_nonbonded_qp_forces_host_v2();
        calc_nonbonded_ww_forces_host_v2();
        calc_nonbonded_qw_forces_host_v2();
    }

    host.E_bond_p.Ubond = calc_bond_forces_host(0, host.n_bonds_solute);
    host.E_bond_w.Ubond = calc_bond_forces_host(host.n_bonds_solute, host.n_bonds);

    host.E_bond_p.Uangle = calc_angle_forces_host(0, host.n_angles_solute);
    host.E_bond_w.Uangle = calc_angle_forces_host(host.n_angles_solute, host.n_angles);

    host.E_bond_p.Utor = calc_torsion_forces_host(0, host.n_torsions_solute);
    host.E_bond_w.Utor = calc_torsion_forces_host(host.n_torsions_solute, host.n_torsions);

    host.E_bond_p.Uimp = calc_improper2_forces_host(0, host.n_impropers_solute);
    host.E_bond_w.Uimp = calc_improper2_forces_host(host.n_impropers_solute, host.n_impropers);

    calc_pshell_forces_host();
    calc_restrseq_forces_host();
    calc_restrpos_forces_host();
    calc_restrdis_forces_host();
    calc_restrang_force_host();
    calc_restrwall_forces_host();

    if (host.n_waters > 0) {
        calc_radix_water_forces_host();
        if (host.md.polarisation) {
            calc_polx_water_forces_host(iteration);
        }
    }

    calc_q_q_nonbonded_forces_host_fallback(host);
    calc_nonbonded_14_forces_host();
    calc_q_bonded_forces_host_fallback(host);
}

void CudaHandler::calc_nonbonded_forces() {
}

void CudaHandler::calc_temperature() {
    ::calc_temperature();
}

void CudaHandler::calc_leapfrog() {
    calc_leapfrog_host();
}

void CudaHandler::prepare_force_dump() {
    auto& host = Context::instance();
    std::vector<fixed_dvel_t> fixed_force_buffer(host.n_atoms);
    host.fixed_dvelocities->download(fixed_force_buffer.data());
    auto* fixed_forces = fixed_force_buffer.data();
    auto* forces = host.dvelocities->cpu_data_p;
    for (int i = 0; i < host.n_atoms; i++) {
        forces[i].x = fixed_to_force(force_fixed_from_storage(fixed_forces[i].x));
        forces[i].y = fixed_to_force(force_fixed_from_storage(fixed_forces[i].y));
        forces[i].z = fixed_to_force(force_fixed_from_storage(fixed_forces[i].z));
    }
}

void CudaHandler::reset_energies() {
    Handler::reset_energies();
    Context::instance().cuda_reset_energies();
}
