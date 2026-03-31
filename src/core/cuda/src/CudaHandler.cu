#include <iostream>

#include "cpu/include/context.h"
#include "cuda/include/CudaHandler.cuh"
#include "cuda/include/CudaAngleForce.cuh"
#include "cuda/include/CudaBondForce.cuh"
#include "cuda/include/CudaContext.cuh"
#include "cuda/include/CudaImproper2Force.cuh"
#include "cuda/include/CudaLeapfrog.cuh"
#include "cuda/include/CudaNonbondedForce.cuh"
#include "cuda/include/CudaNonbondedPPForce.cuh"
#include "cuda/include/CudaNonbondedPWForce.cuh"
#include "cuda/include/CudaNonbondedQPForce.cuh"
#include "cuda/include/CudaNonbondedQQForce.cuh"
#include "cuda/include/CudaNonbondedQWForce.cuh"
#include "cuda/include/CudaNonbondedWWForce.cuh"
#include "cuda/include/CudaPolxWaterForce.cuh"
#include "cuda/include/CudaPshellForce.cuh"
#include "cuda/include/CudaRadixWaterForce.cuh"
#include "cuda/include/CudaRestrangForce.cuh"
#include "cuda/include/CudaRestrdisForce.cuh"
#include "cuda/include/CudaRestrposForce.cuh"
#include "cuda/include/CudaRestrseqForce.cuh"
#include "cuda/include/CudaRestrwallForce.cuh"
#include "cuda/include/CudaShakeConstraints.cuh"
#include "cuda/include/CudaTemperature.cuh"
#include "cuda/include/CudaTorsionForce.cuh"
#include "qatoms.h"

void CudaHandler::initialize() {
    if (!initialized_) {
        ctx_.init();
        init_angle_force_kernel_data();
        init_bond_force_kernel_data();
        init_improper2_force_kernel_data();
        init_leapfrog_kernel_data();
        init_nonbonded_force_kernel_data();
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
        init_shake_constraints_kernel_data();
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

void CudaHandler::run(int num_iterations) {
    for (int i = 0; i < num_iterations; ++i) {
        run_iteration(i);
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

    // todo: No Q-atom bonded iteration now...
}

void CudaHandler::calc_nonbonded_forces() {
    calc_nonbonded_qp_forces_host_v2();
    calc_nonbonded_pp_forces_host_v2();
    calc_nonbonded_ww_forces_host_v2();
    calc_nonbonded_pw_forces_host_v2();
    calc_nonbonded_qw_forces_host_v2();
    calc_nonbonded_qq_forces_host();
}


void CudaHandler::reset_energies() {
    auto& host = Context::instance();
    host.E_total.Upot = 0;
    host.E_bond_p.Uangle = 0;
    host.E_bond_p.Ubond = 0;
    host.E_bond_p.Utor = 0;
    host.E_bond_p.Uimp = 0;
    host.E_bond_w.Uangle = 0;
    host.E_bond_w.Ubond = 0;
    host.E_bond_w.Utor = 0;
    host.E_bond_w.Uimp = 0;
    host.E_bond_q.Uangle = 0;
    host.E_bond_q.Ubond = 0;
    host.E_bond_q.Utor = 0;
    host.E_bond_q.Uimp = 0;
    host.E_nonbond_pp.Ucoul = 0;
    host.E_nonbond_pp.Uvdw = 0;
    host.E_nonbond_pw.Ucoul = 0;
    host.E_nonbond_pw.Uvdw = 0;
    host.E_nonbond_ww.Ucoul = 0;
    host.E_nonbond_ww.Uvdw = 0;
    host.E_nonbond_qx.Ucoul = 0;
    host.E_nonbond_qx.Uvdw = 0;
    host.E_restraint.Uradx = 0;
    host.E_restraint.Upolx = 0;
    host.E_restraint.Ufix = 0;
    host.E_restraint.Ushell = 0;
    host.E_restraint.Upres = 0;
    host.E_restraint.Urestr = 0;
    for (auto& dvel : host.dvelocities) {
        dvel.x = 0;
        dvel.y = 0;
        dvel.z = 0;
    }
    for (int state = 0; state < host.n_lambdas; state++) {
        host.EQ_bond[state].Uangle = 0;
        host.EQ_bond[state].Ubond = 0;
        host.EQ_bond[state].Utor = 0;
        host.EQ_bond[state].Uimp = 0;
        host.EQ_nonbond_qq[state].Ucoul = 0;
        host.EQ_nonbond_qq[state].Uvdw = 0;
        host.EQ_nonbond_qp[state].Ucoul = 0;
        host.EQ_nonbond_qp[state].Uvdw = 0;
        host.EQ_nonbond_qw[state].Ucoul = 0;
        host.EQ_nonbond_qw[state].Uvdw = 0;
        host.EQ_nonbond_qx[state].Ucoul = 0;
        host.EQ_nonbond_qx[state].Uvdw = 0;
        host.EQ_total[state].Utot = 0;
        host.EQ_restraint[state].Urestr = 0;
    }
}

void CudaHandler::run_iteration(int iteration) {
    printf("================================================\n");
    if (iteration > 0) {
        printf("== STEP %d\n", iteration);
    } else {
        printf("== INITIAL ENERGIES\n");
    }
    printf("================================================\n");
    reset_energies();
    // Reset host and device accumulators
    CudaContext& ctx = CudaContext::instance();
    ctx.reset_energies();

    // 1. temperature calculation
    calc_temperature_host();

    // 2. bonded forces and some constraints calculations
    calc_internal_forces(iteration);

    // 3. nonbonded forces
    calc_nonbonded_forces();

    // 5. leapfrog integration
    calc_leapfrog_host();

    // 6. update temperature after integration.
    // todo: Why calculate temperature again here?
    calc_temperature_host();

    // 7. update totals and print outputs like CPU path
    update_energy_totals();
    print_outputs(iteration);
    // Append output files
    write_coords(iteration);
    write_velocities(iteration);
    write_energies(iteration);    
}

void CudaHandler::update_energy_totals() {
    auto& host = Context::instance();
    for (int state = 0; state < host.n_lambdas; state++) {
        if (host.lambdas[state] == 0) {
            host.EQ_bond[state].Uangle = 0;
            host.EQ_bond[state].Ubond = 0;
            host.EQ_bond[state].Utor = 0;
            host.EQ_bond[state].Uimp = 0;
            host.EQ_nonbond_qq[state].Ucoul = 0;
            host.EQ_nonbond_qq[state].Uvdw = 0;
            host.EQ_nonbond_qp[state].Ucoul = 0;
            host.EQ_nonbond_qp[state].Uvdw = 0;
            host.EQ_nonbond_qw[state].Ucoul = 0;
            host.EQ_nonbond_qw[state].Uvdw = 0;
            host.EQ_restraint[state].Urestr = 0;
        }

        host.EQ_nonbond_qx[state].Ucoul = host.EQ_nonbond_qq[state].Ucoul + host.EQ_nonbond_qp[state].Ucoul + host.EQ_nonbond_qw[state].Ucoul;
        host.EQ_nonbond_qx[state].Uvdw = host.EQ_nonbond_qq[state].Uvdw + host.EQ_nonbond_qp[state].Uvdw + host.EQ_nonbond_qw[state].Uvdw;

        host.EQ_total[state].Utot = host.EQ_bond[state].Ubond + host.EQ_bond[state].Uangle + host.EQ_bond[state].Utor + host.EQ_bond[state].Uimp +
                                    host.EQ_nonbond_qx[state].Ucoul + host.EQ_nonbond_qx[state].Uvdw + host.EQ_restraint[state].Urestr;

        host.E_bond_q.Ubond += host.EQ_bond[state].Ubond * host.lambdas[state];
        host.E_bond_q.Uangle += host.EQ_bond[state].Uangle * host.lambdas[state];
        host.E_bond_q.Utor += host.EQ_bond[state].Utor * host.lambdas[state];
        host.E_bond_q.Uimp += host.EQ_bond[state].Uimp * host.lambdas[state];
        host.E_nonbond_qx.Ucoul += host.EQ_nonbond_qx[state].Ucoul * host.lambdas[state];
        host.E_nonbond_qx.Uvdw += host.EQ_nonbond_qx[state].Uvdw * host.lambdas[state];

        // Update protein restraint energies with an average of all states
        host.E_restraint.Upres += host.EQ_restraint[state].Urestr * host.lambdas[state];
    }

    // Update totals
    host.E_restraint.Urestr = host.E_restraint.Uradx + host.E_restraint.Upolx + host.E_restraint.Ushell + host.E_restraint.Ufix + host.E_restraint.Upres;
    host.E_total.Upot = host.E_bond_p.Ubond + host.E_bond_w.Ubond + host.E_bond_p.Uangle + host.E_bond_w.Uangle + host.E_bond_p.Utor + host.E_bond_p.Uimp +
                        host.E_nonbond_pp.Ucoul + host.E_nonbond_pp.Uvdw + host.E_nonbond_pw.Ucoul + host.E_nonbond_pw.Uvdw + host.E_nonbond_ww.Ucoul +
                        host.E_nonbond_ww.Uvdw + host.E_bond_q.Ubond + host.E_bond_q.Uangle + host.E_bond_q.Utor + host.E_bond_q.Uimp +
                        host.E_nonbond_qx.Ucoul + host.E_nonbond_qx.Uvdw + host.E_restraint.Urestr;
    host.E_total.Utot = host.E_total.Upot + host.E_total.Ukin;
}

void CudaHandler::print_outputs(int iteration) {
    auto& host = Context::instance();
    printf("[temperature]\n");
    printf("Temp\t%f\n", host.Temp);
    printf("\n");

    printf("[bonded]\n");
    printf("type\tUbond\tUangle\tUtor\tUimp\n");
    printf("p\t%f\t%f\t%f\t%f\n", host.E_bond_p.Ubond, host.E_bond_p.Uangle, host.E_bond_p.Utor, host.E_bond_p.Uimp);
    printf("w\t%f\t%f\t%f\t%f\n", host.E_bond_w.Ubond, host.E_bond_w.Uangle, host.E_bond_w.Utor, host.E_bond_w.Uimp);
    printf("qp\t%f\t%f\t%f\t%f\n", host.E_bond_q.Ubond, host.E_bond_q.Uangle, host.E_bond_q.Utor, host.E_bond_q.Uimp);
    printf("\n");

    printf("[nonbonded]\n");
    printf("type\tUcoul\tUvdw\n");
    printf("pp\t%f\t%f\n", host.E_nonbond_pp.Ucoul, host.E_nonbond_pp.Uvdw);
    printf("pw\t%f\t%f\n", host.E_nonbond_pw.Ucoul, host.E_nonbond_pw.Uvdw);
    printf("ww\t%f\t%f\n", host.E_nonbond_ww.Ucoul, host.E_nonbond_ww.Uvdw);
    printf("qx\t%f\t%f\n", host.E_nonbond_qx.Ucoul, host.E_nonbond_qx.Uvdw);
    printf("\n");

    printf("[restraint]\n");
    printf("Uradx\t%f\n", host.E_restraint.Uradx);
    printf("Upolx\t%f\n", host.E_restraint.Upolx);
    printf("Ushell\t%f\n", host.E_restraint.Ushell);
    printf("Ufix\t%f\n", host.E_restraint.Ufix);
    printf("Upres\t%f\n", host.E_restraint.Upres);
    printf("Total\t%f\n", host.E_restraint.Urestr);
    printf("\n");

    printf("[q-energies]\n");
    printf("lambda\tSUM\tUbond\tUangle\tUtor\tUimp\tUcoul\tUvdw\tUrestr\n");
    for (int state = 0; state < host.n_lambdas; state++) {
        printf("%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\n", host.lambdas[state], host.EQ_total[state].Utot, host.EQ_bond[state].Ubond,
               host.EQ_bond[state].Uangle, host.EQ_bond[state].Utor, host.EQ_bond[state].Uimp, host.EQ_nonbond_qx[state].Ucoul, host.EQ_nonbond_qx[state].Uvdw, host.EQ_restraint[state].Urestr);
    }
    printf("\n");

    printf("[total]\n");
    printf("Ukin\t%f\n", host.E_total.Ukin);
    printf("Upot\t%f\n", host.E_total.Upot);
    printf("Utot\t%f\n", host.E_total.Utot);
    printf("\n");
}
