#include "water_boundary_force.h"

#include <math.h>

#include "constants.h"
#include "q_math.h"

void WaterBoundaryForce::init(Context& ctx, const Temperature& temperature, const ParseResult& parsed) {
    temperature_ = &temperature;
    if (ctx.n_waters() > 0) {
        const double crgQtot = build_water_sphere(ctx);
        build_wshells(ctx, crgQtot, parsed);
    }
    init_backend(ctx);
}

double WaterBoundaryForce::build_water_sphere(Context& ctx) {
    double crgQtot = 0;
    if (ctx.md.charge_correction && ctx.n_lambdas() > 0) {
        auto* lambdas = ctx.lambdas->cpu_data_p;
        for (int state = 0; state < ctx.n_lambdas(); state++) {
            for (int qi = 0; qi < ctx.n_qatoms(); qi++) {
                crgQtot += ctx.q_charges[qi + ctx.n_qatoms() * state].charge * lambdas[state];
            }
        }
    }
    /*
    U = Dwmz·(e^(2·awmz·db) − 2·e^(awmz·db)),
    a well of depth Dwmz at the wall, and awmz controls how fast that surface attraction decays as you move inward.
    (awmz has units of 1/length because awmz·db must be dimensionless with db in Å.)
    */
    data_.Dwmz = 0.26 * exp(-0.19 * (ctx.topo.solvent_radius - 15)) + 0.74;
    data_.awmz = 0.2 / (1 + exp(0.4 * (ctx.topo.solvent_radius - 25))) + 0.3;

    printf("Dwmz = %f, awmz = %f\n", data_.Dwmz, data_.awmz);
    return crgQtot;
}

void WaterBoundaryForce::build_wshells(Context& ctx, double crgQtot, const ParseResult& parsed) {
    bool run_gpu = ctx.command_info.requested_gpu;

    const auto& restart_theta_corr = parsed.restart_theta_corr;
    const auto& water_bond = parsed.bonds[parsed.n_bonds_solute];
    const auto& water_angle = parsed.angles[parsed.n_angles_solute];
    const auto& water_bond_parameter = parsed.cbonds[water_bond.code - 1];
    const auto& water_angle_parameter = parsed.cangles[water_angle.code - 1];
    const int oxygen = ctx.n_atoms_solute;
    const charge_t& oxygen_charge = parsed.charges[oxygen];
    const ccharge_t& oxygen_charge_parameter = parsed.ccharges[oxygen_charge.code - 1];
    double crg_ow = oxygen_charge_parameter.charge;
    const double mu_w = -crg_ow * water_bond_parameter.b0 * cos(deg2rad(water_angle_parameter.th0) / 2);

    const double drs = wpolr_layer / drouter;
    data_.n_shells = (int)floor(-0.5 + sqrt(2 * drs + 0.25));
    data_.wshells = std::make_unique<HostDeviceBuffer<shell_t>>(data_.n_shells, true, run_gpu);

    auto* wshells = data_.wshells->cpu_data_p;
    printf("n_shells = %d\n", data_.n_shells);

    double router = ctx.topo.solvent_radius;
    data_.n_max_inshell = 0;
    for (int i = 0; i < data_.n_shells; i++) {
        wshells[i].n_inshell = 0;
        wshells[i].theta_corr = 0;
        wshells[i].avtheta = 0;
        wshells[i].avn_inshell = 0;
        wshells[i].router = router;
        const double dr = drouter * (i + 1);
        const double ri = router - dr;

        wshells[i].dr = dr;
        const double vshell = pow(router, 3) - pow(ri, 3);
        const int n_inshell = static_cast<int>(floor(4 * M_PI / 3 * vshell * rho_water));
        if (n_inshell > data_.n_max_inshell) {
            data_.n_max_inshell = n_inshell;
        }

        const double rshell = pow(0.5 * (pow(router, 3) + pow(ri, 3)), 1.0 / 3.0);
        wshells[i].cstb = crgQtot * 0.98750 / (rho_water * mu_w * 4 * M_PI * pow(rshell, 2));  // 0.98750 = (1-1/epsilon) for water
        router -= dr;
    }

    if (restart_theta_corr.size() == static_cast<size_t>(data_.n_shells)) {
        for (int i = 0; i < data_.n_shells; i++) {
            wshells[i].theta_corr = restart_theta_corr[i];
        }
    }

    for (int i = 0; i < data_.n_shells; i++) {
        printf("shell %d: (%f, %f)\n",  // verbatim (see note)
               i, wshells[i].router, wshells[i].router - wshells[i].dr);
    }
    // data_.n_max_inshell *= 1.5;  // take largest and add some extra
    data_.n_max_inshell = ctx.n_waters();
    data_.theta = std::make_unique<HostDeviceBuffer<double>>(ctx.n_waters(), true, run_gpu);
    data_.list_sh = std::make_unique<HostDeviceBuffer<int>>(data_.n_max_inshell * data_.n_shells, true, run_gpu);

    data_.theta->zero();
    data_.list_sh->zero();

    if (run_gpu) data_.wshells->upload();
}

void WaterBoundaryForce::sync_for_output(
    Context& ctx) {
    if (!data_.wshells ||
        data_.n_shells <= 0) {
        ctx.wshells.clear();
        return;
    }

    const shell_t* source =
        data_.wshells->cpu_data_p;

    ctx.wshells.assign(
        source,
        source + data_.n_shells);
}