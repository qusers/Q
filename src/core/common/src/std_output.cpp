#include "std_output.h"

#include <cstdio>

void StdOutput::output_trajectory(Context& ctx, int iteration) {
}

void StdOutput::output_energy(Context& ctx, int iteration) {
     if (iteration > 0 && (ctx.md.output <= 0 || iteration % ctx.md.output != 0)) return;
    auto* lambdas = ctx.lambdas ? ctx.lambdas->cpu_data_p : nullptr;
    auto& energy = ctx.energy.data();

    std::printf("================================================\n");
    if (iteration > 0) {
        std::printf("== STEP %d\n", iteration);
    } else if (iteration == 0) {
        std::printf("== INITIAL ENERGIES\n");
    } else {
        std::printf("== FINAL ENERGIES summary\n");
    }
    std::printf("================================================\n");

    std::printf("[temperature]\n");
    std::printf("Temp\t%f\n", ctx.Temp);
    std::printf("\n");

    std::printf("[bonded]\n");
    std::printf("type\tUbond\tUangle\tUtor\tUimp\n");
    std::printf("p\t%f\t%f\t%f\t%f\n", energy.bond_p.Ubond, energy.bond_p.Uangle, energy.bond_p.Utor, energy.bond_p.Uimp);
    std::printf("w\t%f\t%f\t%f\t%f\n", energy.bond_w.Ubond, energy.bond_w.Uangle, energy.bond_w.Utor, energy.bond_w.Uimp);
    std::printf("qp\t%f\t%f\t%f\t%f\n", energy.bond_q.Ubond, energy.bond_q.Uangle, energy.bond_q.Utor, energy.bond_q.Uimp);
    std::printf("\n");

    std::printf("[nonbonded]\n");
    std::printf("type\tUcoul\tUvdw\n");
    std::printf("pp\t%f\t%f\n", energy.nb_pp.Ucoul, energy.nb_pp.Uvdw);
    std::printf("pw\t%f\t%f\n", energy.nb_pw.Ucoul, energy.nb_pw.Uvdw);
    std::printf("ww\t%f\t%f\n", energy.nb_ww.Ucoul, energy.nb_ww.Uvdw);
    std::printf("qx\t%f\t%f\n", energy.nb_qx.Ucoul, energy.nb_qx.Uvdw);
    std::printf("\n");

    std::printf("[restraint]\n");
    std::printf("Uradx\t%f\n", energy.restraint.Uradx);
    std::printf("Upolx\t%f\n", energy.restraint.Upolx);
    std::printf("Ushell\t%f\n", energy.restraint.Ushell);
    std::printf("Ufix\t%f\n", energy.restraint.Ufix);
    std::printf("Upres\t%f\n", energy.restraint.Upres);
    std::printf("Total\t%f\n", energy.restraint.Urestr);
    std::printf("\n");

    std::printf("[q-energies]\n");
    std::printf("lambda\tSUM\tUbond\tUangle\tUtor\tUimp\tUcoul\tUvdw\tUrestr\n");
    for (int state = 0; state < ctx.n_lambdas; state++) {
        std::printf("%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\n",
                    lambdas[state],
                    energy.eq_total[state],
                    energy.eq_bond[state].Ubond,
                    energy.eq_bond[state].Uangle,
                    energy.eq_bond[state].Utor,
                    energy.eq_bond[state].Uimp,
                    energy.eq_qx[state].Ucoul,
                    energy.eq_qx[state].Uvdw,
                    energy.eq_restr[state]);
    }
    std::printf("\n");

    std::printf("[total]\n");
    std::printf("Ukin\t%f\n", energy.Ukin);
    std::printf("Upot\t%f\n", energy.Upot);
    std::printf("Utot\t%f\n", energy.Utot);
    std::printf("\n");
}

void StdOutput::finish(Context& ctx) {
    output_energy(ctx, -1);
}