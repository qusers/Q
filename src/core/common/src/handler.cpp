#include "handler.h"

#include "output.h"

void Handler::run_iteration(int iteration) {
    printf("================================================\n");
    if (iteration > 0) {
        printf("== STEP %d\n", iteration);
    } else {
        printf("== INITIAL ENERGIES\n");
    }
    printf("================================================\n");

    reset_energies();

    // 1. temperature calculation
    calc_temperature();

    // 2. bonded forces and some constraints calculations
    calc_internal_forces(iteration);

    // 3. nonbonded forces
    calc_nonbonded_forces();

    // 5. leapfrog integration
    calc_leapfrog();

    // 6. update temperature after integration.
    // todo: Why calculate temperature again here?
    calc_temperature();

    // 7. update totals energy
    update_energy_totals();

    // 8. print output to files
    print_outputs(iteration);
}

void Handler::run(int num_iterations) {
    for (int i = 0; i < num_iterations; i++) {
        run_iteration(i);
    }
}

void Handler::update_energy_totals() {
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

void Handler::print_outputs(int iteration) {
    print_energies();
    write_coords(iteration);
    write_velocities(iteration);
    write_energies(iteration);
}

void Handler::reset_energies() {
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
