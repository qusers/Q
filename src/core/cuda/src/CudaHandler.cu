#include <iostream>

#include "CudaHandler.cuh"
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
    E_bond_p.Uangle = calc_angle_forces_host(0, n_angles_solute);
    E_bond_w.Uangle = calc_angle_forces_host(n_angles_solute, n_angles);

    E_bond_p.Ubond = calc_bond_forces_host(0, n_bonds_solute);
    E_bond_w.Ubond = calc_bond_forces_host(n_bonds_solute, n_bonds);

    E_bond_p.Utor = calc_torsion_forces_host(0, n_torsions_solute);
    E_bond_w.Utor = calc_torsion_forces_host(n_torsions_solute, n_torsions);

    E_bond_p.Uimp = calc_improper2_forces_host(0, n_impropers_solute);
    E_bond_w.Uimp = calc_improper2_forces_host(n_impropers_solute, n_impropers);

    if (n_waters > 0) {
        calc_radix_water_forces_host();
        if (md.polarisation) {
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
    E_total.Upot = 0;
    E_bond_p.Uangle = 0;
    E_bond_p.Ubond = 0;
    E_bond_p.Utor = 0;
    E_bond_p.Uimp = 0;
    E_bond_w.Uangle = 0;
    E_bond_w.Ubond = 0;
    E_bond_w.Utor = 0;
    E_bond_w.Uimp = 0;
    E_bond_q.Uangle = 0;
    E_bond_q.Ubond = 0;
    E_bond_q.Utor = 0;
    E_bond_q.Uimp = 0;
    E_nonbond_pp.Ucoul = 0;
    E_nonbond_pp.Uvdw = 0;
    E_nonbond_pw.Ucoul = 0;
    E_nonbond_pw.Uvdw = 0;
    E_nonbond_ww.Ucoul = 0;
    E_nonbond_ww.Uvdw = 0;
    E_nonbond_qx.Ucoul = 0;
    E_nonbond_qx.Uvdw = 0;
    E_restraint.Uradx = 0;
    E_restraint.Upolx = 0;
    E_restraint.Ufix = 0;
    E_restraint.Ushell = 0;
    E_restraint.Upres = 0;
    E_restraint.Urestr = 0;
    for (int state = 0; state < n_lambdas; state++) {
        EQ_bond[state].Uangle = 0;
        EQ_bond[state].Ubond = 0;
        EQ_bond[state].Utor = 0;
        EQ_bond[state].Uimp = 0;
        EQ_nonbond_qq[state].Ucoul = 0;
        EQ_nonbond_qq[state].Uvdw = 0;
        EQ_nonbond_qp[state].Ucoul = 0;
        EQ_nonbond_qp[state].Uvdw = 0;
        EQ_nonbond_qw[state].Ucoul = 0;
        EQ_nonbond_qw[state].Uvdw = 0;
        EQ_restraint[state].Urestr = 0;
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
}

void CudaHandler::update_energy_totals() {
    for (int state = 0; state < n_lambdas; state++) {

        EQ_nonbond_qx[state].Ucoul = EQ_nonbond_qq[state].Ucoul + EQ_nonbond_qp[state].Ucoul + EQ_nonbond_qw[state].Ucoul;
        EQ_nonbond_qx[state].Uvdw = EQ_nonbond_qq[state].Uvdw + EQ_nonbond_qp[state].Uvdw + EQ_nonbond_qw[state].Uvdw;

        EQ_total[state].Utot = EQ_bond[state].Ubond + EQ_bond[state].Uangle + EQ_bond[state].Utor + EQ_bond[state].Uimp + EQ_nonbond_qx[state].Ucoul + EQ_nonbond_qx[state].Uvdw + EQ_restraint[state].Urestr;

        E_bond_q.Ubond += EQ_bond[state].Ubond * lambdas[state];
        E_bond_q.Uangle += EQ_bond[state].Uangle * lambdas[state];
        E_bond_q.Utor += EQ_bond[state].Utor * lambdas[state];
        E_bond_q.Uimp += EQ_bond[state].Uimp * lambdas[state];
        E_nonbond_qx.Ucoul += EQ_nonbond_qx[state].Ucoul * lambdas[state];
        E_nonbond_qx.Uvdw += EQ_nonbond_qx[state].Uvdw * lambdas[state];

        // Update protein restraint energies with an average of all states
        E_restraint.Upres += EQ_restraint[state].Urestr * lambdas[state];
    }

    // Update totals
    E_restraint.Urestr = E_restraint.Uradx + E_restraint.Upolx + E_restraint.Ushell + E_restraint.Ufix + E_restraint.Upres;
    E_total.Upot = E_bond_p.Ubond + E_bond_w.Ubond + E_bond_p.Uangle + E_bond_w.Uangle + E_bond_p.Utor + E_bond_p.Uimp + E_nonbond_pp.Ucoul + E_nonbond_pp.Uvdw + E_nonbond_pw.Ucoul + E_nonbond_pw.Uvdw + E_nonbond_ww.Ucoul + E_nonbond_ww.Uvdw + E_bond_q.Ubond + E_bond_q.Uangle + E_bond_q.Utor + E_bond_q.Uimp + E_nonbond_qx.Ucoul + E_nonbond_qx.Uvdw + E_restraint.Urestr;
    E_total.Utot = E_total.Upot + E_total.Ukin;
}

void CudaHandler::print_outputs(int iteration) {
    printf("[temperature]\n");
    printf("Temp\t%f\n", Temp);
    printf("\n");

    printf("[bonded]\n");
    printf("type\tUbond\tUangle\tUtor\tUimp\n");
    printf("p\t%f\t%f\t%f\t%f\n", E_bond_p.Ubond, E_bond_p.Uangle, E_bond_p.Utor, E_bond_p.Uimp);
    printf("w\t%f\t%f\t%f\t%f\n", E_bond_w.Ubond, E_bond_w.Uangle, E_bond_w.Utor, E_bond_w.Uimp);
    printf("qp\t%f\t%f\t%f\t%f\n", E_bond_q.Ubond, E_bond_q.Uangle, E_bond_q.Utor, E_bond_q.Uimp);
    printf("\n");

    printf("[nonbonded]\n");
    printf("type\tUcoul\tUvdw\n");
    printf("pp\t%f\t%f\n", E_nonbond_pp.Ucoul, E_nonbond_pp.Uvdw);
    printf("pw\t%f\t%f\n", E_nonbond_pw.Ucoul, E_nonbond_pw.Uvdw);
    printf("ww\t%f\t%f\n", E_nonbond_ww.Ucoul, E_nonbond_ww.Uvdw);
    printf("qx\t%f\t%f\n", E_nonbond_qx.Ucoul, E_nonbond_qx.Uvdw);
    printf("\n");

    printf("[restraint]\n");
    printf("Uradx\t%f\n", E_restraint.Uradx);
    printf("Upolx\t%f\n", E_restraint.Upolx);
    printf("Ushell\t%f\n", E_restraint.Ushell);
    printf("Ufix\t%f\n", E_restraint.Ufix);
    printf("Upres\t%f\n", E_restraint.Upres);
    printf("Total\t%f\n", E_restraint.Urestr);
    printf("\n");

    printf("[q-energies]\n");
    printf("lambda\tSUM\tUbond\tUangle\tUtor\tUimp\tUcoul\tUvdw\tUrestr\n");
    for (int state = 0; state < n_lambdas; state++) {
        printf("%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\n", lambdas[state], EQ_total[state].Utot, EQ_bond[state].Ubond,
               EQ_bond[state].Uangle, EQ_bond[state].Utor, EQ_bond[state].Uimp, EQ_nonbond_qx[state].Ucoul, EQ_nonbond_qx[state].Uvdw, EQ_restraint[state].Urestr);
    }
    printf("\n");

    printf("[total]\n");
    printf("Ukin\t%f\n", E_total.Ukin);
    printf("Upot\t%f\n", E_total.Upot);
    printf("Utot\t%f\n", E_total.Utot);
    printf("\n");
}
