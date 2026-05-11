#include "std_output.h"

#include <cstdio>

void StdOutput::output_trajectory(Context& ctx, int iteration) {
    (void) ctx;
    (void) iteration;
}

void StdOutput::output_energy(Context& ctx, int iteration) {
    (void) iteration;

    auto* lambdas = ctx.lambdas ? ctx.lambdas->cpu_data_p : nullptr;
    auto* EQ_restraint = ctx.EQ_restraint->cpu_data_p;

    std::printf("[temperature]\n");
    std::printf("Temp\t%f\n", ctx.Temp);
    std::printf("\n");

    std::printf("[bonded]\n");
    std::printf("type\tUbond\tUangle\tUtor\tUimp\n");
    std::printf("p\t%f\t%f\t%f\t%f\n", ctx.E_bond_p.Ubond, ctx.E_bond_p.Uangle, ctx.E_bond_p.Utor, ctx.E_bond_p.Uimp);
    std::printf("w\t%f\t%f\t%f\t%f\n", ctx.E_bond_w.Ubond, ctx.E_bond_w.Uangle, ctx.E_bond_w.Utor, ctx.E_bond_w.Uimp);
    std::printf("qp\t%f\t%f\t%f\t%f\n", ctx.E_bond_q.Ubond, ctx.E_bond_q.Uangle, ctx.E_bond_q.Utor, ctx.E_bond_q.Uimp);
    std::printf("\n");

    std::printf("[nonbonded]\n");
    std::printf("type\tUcoul\tUvdw\n");
    std::printf("pp\t%f\t%f\n", ctx.E_nonbond_pp.Ucoul, ctx.E_nonbond_pp.Uvdw);
    std::printf("pw\t%f\t%f\n", ctx.E_nonbond_pw.Ucoul, ctx.E_nonbond_pw.Uvdw);
    std::printf("ww\t%f\t%f\n", ctx.E_nonbond_ww.Ucoul, ctx.E_nonbond_ww.Uvdw);
    std::printf("qx\t%f\t%f\n", ctx.E_nonbond_qx.Ucoul, ctx.E_nonbond_qx.Uvdw);
    std::printf("\n");

    std::printf("[restraint]\n");
    std::printf("Uradx\t%f\n", ctx.E_restraint.Uradx);
    std::printf("Upolx\t%f\n", ctx.E_restraint.Upolx);
    std::printf("Ushell\t%f\n", ctx.E_restraint.Ushell);
    std::printf("Ufix\t%f\n", ctx.E_restraint.Ufix);
    std::printf("Upres\t%f\n", ctx.E_restraint.Upres);
    std::printf("Total\t%f\n", ctx.E_restraint.Urestr);
    std::printf("\n");

    std::printf("[q-energies]\n");
    std::printf("lambda\tSUM\tUbond\tUangle\tUtor\tUimp\tUcoul\tUvdw\tUrestr\n");
    for (int state = 0; state < ctx.n_lambdas; state++) {
        std::printf("%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\n",
                    lambdas[state],
                    ctx.EQ_total[state].Utot,
                    ctx.EQ_bond[state].Ubond,
                    ctx.EQ_bond[state].Uangle,
                    ctx.EQ_bond[state].Utor,
                    ctx.EQ_bond[state].Uimp,
                    ctx.EQ_nonbond_qx[state].Ucoul,
                    ctx.EQ_nonbond_qx[state].Uvdw,
                    EQ_restraint[state].Urestr);
    }
    std::printf("\n");

    std::printf("[total]\n");
    std::printf("Ukin\t%f\n", ctx.E_total.Ukin);
    std::printf("Upot\t%f\n", ctx.E_total.Upot);
    std::printf("Utot\t%f\n", ctx.E_total.Utot);
    std::printf("\n");
}
