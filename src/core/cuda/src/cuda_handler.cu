#include <iostream>

#include "common/include/context.h"
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
#include "cuda/include/cuda_shake_constraints.cuh"
#include "cuda/include/cuda_temperature.cuh"
#include "cuda/include/cuda_torsion_force.cuh"

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
        init_shake_constraints_kernel_data(Context::instance());
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
        cleanup_shake_constraints();
        cleanup_temperature();
        cleanup_torsion_force();
        initialized_ = false;
    }
}

void CudaHandler::calc_internal_forces(int iteration) {
    auto& host = Context::instance();

    host.E_bond_p.Uangle = calc_angle_forces_host(0, host.n_angles_solute);
    host.E_bond_w.Uangle = calc_angle_forces_host(host.n_angles_solute, host.n_angles);

    host.E_bond_p.Ubond = calc_bond_forces_host(0, host.n_bonds_solute);
    host.E_bond_w.Ubond = calc_bond_forces_host(host.n_bonds_solute, host.n_bonds);

    host.E_bond_p.Utor = calc_torsion_forces_host(0, host.n_torsions_solute);
    host.E_bond_w.Utor = calc_torsion_forces_host(host.n_torsions_solute, host.n_torsions);

    host.E_bond_p.Uimp = calc_improper2_forces_host(0, host.n_impropers_solute);
    host.E_bond_w.Uimp = calc_improper2_forces_host(host.n_impropers_solute, host.n_impropers);

    if (host.n_waters > 0) {
        calc_radix_water_forces_host();
        if (host.md.polarisation) {
            calc_polx_water_forces_host(iteration);
        }
    }

    calc_pshell_forces_host();
    calc_restrseq_forces_host();
    calc_restrdis_forces_host();
    calc_restrpos_forces_host();
    calc_restrang_force_host();
    calc_restrwall_forces_host();
}

void CudaHandler::calc_nonbonded_forces() {
    auto& host = Context::instance();
    calc_nonbonded_qp_forces_host_v2();
    calc_nonbonded_pp_forces_host_v2();
    calc_nonbonded_ww_forces_host_v2();
    calc_nonbonded_pw_forces_host_v2();
    calc_nonbonded_qw_forces_host_v2();
    calc_nonbonded_qq_forces_host();
    calc_nonbonded_14_forces_host();
}

void CudaHandler::calc_temperature() {
    ::calc_temperature_host();
}

void CudaHandler::calc_leapfrog() {
    calc_leapfrog_host();
}

void CudaHandler::reset_energies() {
    Handler::reset_energies();
    Context::instance().cuda_reset_energies();
}
