#include "init.h"

#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>
#include <unistd.h>

#include <algorithm>

#include "constants.h"
#include "context.h"
#include "cpu_handler.h"
#include "cpu_utils.h"
#include "parse.h"
#include "shake.h"

// Remove bonds, angles, torsions and impropers which are excluded or changed in the FEP file
void exclude_qatom_definitions() {
    auto& ctx = Context::instance();
    int excluded;
    int ai = 0, bi = 0, ii = 0, ti = 0;
    int qai = 0, qbi = 0, qii = 0, qti = 0;

    excluded = 0;
    int solute_excluded = 0;

    auto& bonds = ctx.bonds;
    auto& angles = ctx.angles;
    auto& torsions = ctx.torsions;

    if (ctx.n_qangles > 0) {
        for (int i = 0; i < ctx.n_angles; i++) {
            if (angles[i].ai == ctx.q_angles[qai].ai && angles[i].aj == ctx.q_angles[qai].aj && angles[i].ak == ctx.q_angles[qai].ak) {
                qai++;
                excluded++;
                if (i < ctx.n_angles_solute) solute_excluded++;
            } else {
                angles[ai] = angles[i];
                ai++;
            }
        }
        ctx.n_angles -= excluded;
        ctx.n_angles_solute -= solute_excluded;
    }

    excluded = 0;
    solute_excluded = 0;
    if (ctx.n_qbonds > 0) {
        for (int i = 0; i < ctx.n_bonds; i++) {
            if (bonds[i].ai == ctx.q_bonds[qbi].ai && bonds[i].aj == ctx.q_bonds[qbi].aj) {
                qbi++;
                excluded++;
                if (i < ctx.n_bonds_solute) solute_excluded++;
            } else {
                bonds[bi] = bonds[i];
                bi++;
            }
        }
        ctx.n_bonds -= excluded;
        ctx.n_bonds_solute -= solute_excluded;
    }

    // excluded = 0;
    // for (int i = 0; i < n_impropers; i++) {
    //     if (impropers[i].ai == q_impropers[qai][0].ai
    //      && impropers[i].aj == q_impropers[qai][0].aj
    //      && impropers[i].ak == q_impropers[qai][0].ak
    //      && impropers[i].al == q_impropers[qai][0].al) {
    //         qii++;
    //         excluded++;
    //     }
    //     else {
    //         impropers[ii] = impropers[i];
    //         ii++;
    //     }
    // }
    // n_impropers -= excluded;

    excluded = 0;
    solute_excluded = 0;
    if (ctx.n_qtorsions > 0) {
        for (int i = 0; i < ctx.n_torsions; i++) {
            if (torsions[i].ai == ctx.q_torsions[qti].ai && torsions[i].aj == ctx.q_torsions[qti].aj && torsions[i].ak == ctx.q_torsions[qti].ak && torsions[i].al == ctx.q_torsions[qti].al) {
                qti++;
                excluded++;
                if (i < ctx.n_torsions_solute) solute_excluded++;
            } else {
                torsions[ti] = torsions[i];
                ti++;
            }
        }
        ctx.n_torsions -= excluded;
        ctx.n_torsions_solute -= solute_excluded;
    }

    // TODO: add exclusion pairs
}

void exclude_all_atoms_excluded_definitions() {
    auto& ctx = Context::instance();
    auto* excluded = ctx.excluded->cpu_data_p;
    int n_excl;
    int ai = 0, bi = 0, ii = 0, ti = 0;
    auto& impropers = ctx.impropers;
    auto& torsions = ctx.torsions;

    // n_excl = 0;
    // for (int i = 0; i < n_angles; i++) {
    //     if (excluded[angles[i].ai - 1]
    //         && excluded[angles[i].aj - 1]
    //         && excluded[angles[i].ak - 1]) {
    //         n_excl++;
    //     }
    //     else {
    //         angles[ai] = angles[i];
    //         ai++;
    //     }
    // }
    // printf("original: %d. # excluded angles: %d\n", n_angles, n_excl);
    // n_angles -= n_excl;

    // n_excl = 0;
    // for (int i = 0; i < n_bonds; i++) {
    //     if (excluded[bonds[i].ai - 1]
    //         && excluded[bonds[i].aj - 1]) {
    //         n_excl++;
    //     }
    //     else {
    //         bonds[bi] = bonds[i];
    //         bi++;
    //     }
    // }
    // printf("original: %d. # excluded bonds: %d\n", n_bonds, n_excl);
    // n_bonds -= n_excl;

    n_excl = 0;
    for (int i = 0; i < ctx.n_impropers; i++) {
        if (excluded[impropers[i].ai - 1] && excluded[impropers[i].aj - 1] && excluded[impropers[i].ak - 1] && excluded[impropers[i].al - 1]) {
            n_excl++;
        } else {
            impropers[ii] = impropers[i];
            ii++;
        }
    }
    printf("original: %d. # excluded impropers: %d\n", ctx.n_impropers, n_excl);
    ctx.n_impropers -= n_excl;

    n_excl = 0;
    for (int i = 0; i < ctx.n_torsions; i++) {
        if (excluded[torsions[i].ai - 1] && excluded[torsions[i].aj - 1] && excluded[torsions[i].ak - 1] && excluded[torsions[i].al - 1]) {
            n_excl++;
        } else {
            torsions[ti] = torsions[i];
            ti++;
        }
    }
    printf("original: %d. # excluded torsions: %d\n", ctx.n_torsions, n_excl);
    ctx.n_torsions -= n_excl;
}

void init_velocities() {
    auto& ctx = Context::instance();
    if (ctx.velocities) {
        if (ctx.command_info.requested_gpu) {
            ctx.velocities->upload();
        }
        return;
    }

    auto& atypes = ctx.atypes->cpu_data_p;
    auto& catypes = ctx.catypes->cpu_data_p;
    ctx.velocities = std::make_unique<HostDeviceBuffer<vel_t>>(ctx.n_atoms, true, ctx.command_info.requested_gpu);
    auto& velocities = ctx.velocities->cpu_data_p;

    // If not previous value set, use a Maxwell distribution to fill velocities
    double kT = Boltz * ctx.md.initial_temperature;
    double sd, mass;
    for (int i = 0; i < ctx.n_atoms; i++) {
        mass = catypes[atypes[i].code - 1].m;
        sd = sqrt(kT / mass);

        velocities[i].x = gauss(0, sd);
        velocities[i].y = gauss(0, sd);
        velocities[i].z = gauss(0, sd);
    }

    if (ctx.command_info.requested_gpu) {
        ctx.velocities->upload();
    }
}

void init_inv_mass() {
    auto& ctx = Context::instance();
    auto& atypes = ctx.atypes->cpu_data_p;
    auto& catypes = ctx.catypes->cpu_data_p;
    ctx.winv = std::make_unique<HostDeviceBuffer<double>>(ctx.n_atoms, true, ctx.command_info.requested_gpu);
    auto* winv = ctx.winv->cpu_data_p;
    for (int ai = 0; ai < ctx.n_atoms; ai++) {
        winv[ai] = 1 / catypes[atypes[ai].code - 1].m;
    }

    if (ctx.command_info.requested_gpu) {
        ctx.winv->upload();
    }
}

/* =============================================
 * == BOUNDARY RESTRAINTS
 * =============================================
 */

void init_water_sphere() {
    auto& ctx = Context::instance();
    ctx.crgQtot = 0.0;
    if (ctx.md.charge_correction && ctx.n_lambdas > 0) {
        auto* lambdas = ctx.lambdas->cpu_data_p;
        for (int state = 0; state < ctx.n_lambdas; state++) {
            for (int qi = 0; qi < ctx.n_qatoms; qi++) {
                ctx.crgQtot += ctx.q_charges[qi + ctx.n_qatoms * state].charge * lambdas[state];
            }
        }
    }

    ctx.Dwmz = 0.26 * exp(-0.19 * (ctx.topo.solvent_radius - 15)) + 0.74;
    ctx.awmz = 0.2 / (1 + exp(0.4 * (ctx.topo.solvent_radius - 25))) + 0.3;

    printf("Dwmz = %f, awmz = %f\n", ctx.Dwmz, ctx.awmz);
}

// ONLY call if there are actually solvent atoms, or get segfaulted
void init_wshells() {
    auto& ctx = Context::instance();
    int n_inshell;
    double drs, router, ri, dr, Vshell, rshell;
    auto& cbonds = ctx.cbonds;
    auto& cangles = ctx.cangles;
    // Match Fortran wat_shells, which uses the last topology water bond/angle codes.
    cbond_t cbondw = cbonds[ctx.n_cbonds - 1];
    cangle_t canglew = cangles[ctx.n_cangles - 1];
    const double crg_ow = ctx.unified_ccharge(ctx.n_atoms_solute, 0).charge;
    const double mu_w = -crg_ow * cbondw.b0 * cos(to_radians(canglew.th0) / 2);

    drs = wpolr_layer / drouter;

    ctx.n_shells = (int)floor(-0.5 + sqrt(2 * drs + 0.25));
    ctx.wshells = std::make_unique<HostDeviceBuffer<shell_t>>(ctx.n_shells, true, ctx.command_info.requested_gpu);
    auto* wshells = ctx.wshells->cpu_data_p;

    printf("n_shells = %d\n", ctx.n_shells);

    router = ctx.topo.solvent_radius;
    ctx.n_max_inshell = 0;

    for (int i = 0; i < ctx.n_shells; i++) {
        wshells[i].n_inshell = 0;
        wshells[i].theta_corr = 0;
        wshells[i].avtheta = 0;
        wshells[i].avn_inshell = 0;
        wshells[i].router = router;
        dr = drouter * (i + 1);
        ri = router - dr;
        wshells[i].dr = dr;
        Vshell = pow(router, 3) - pow(ri, 3);
        n_inshell = (int)floor(4 * M_PI / 3 * Vshell * rho_water);
        if (n_inshell > ctx.n_max_inshell) {
            ctx.n_max_inshell = n_inshell;
        }
        rshell = pow(0.5 * (pow(router, 3) + pow(ri, 3)), 1.0 / 3.0);

        // --- Note below: 0.98750 = (1-1/epsilon) for water
        wshells[i].cstb = ctx.crgQtot * 0.98750 / (rho_water * mu_w * 4 * M_PI * pow(rshell, 2));

        router -= dr;
    }

    if (ctx.restart_theta_corr.size() == static_cast<size_t>(ctx.n_shells)) {
        for (int i = 0; i < ctx.n_shells; i++) {
            wshells[i].theta_corr = ctx.restart_theta_corr[i];
        }
    }

    // rc > wshells[n_shells-1].router - wshells[n_shells-1].dr
    printf("shell 0: (%f, %f). shell 1: (%f, %f). shell 2: (%f, %f).\n", wshells[0].router, wshells[0].router - wshells[0].dr, wshells[1].router, wshells[1].router - wshells[1].dr, wshells[2].router, wshells[2].router - wshells[2].dr);

    ctx.n_max_inshell = ctx.n_waters;  // Make largest a little bigger just in case

    // Initialize arrays needed for bookkeeping
    ctx.theta.assign(ctx.n_waters, 0.0);
    ctx.theta0.assign(ctx.n_waters, 0.0);
    ctx.tdum.assign(ctx.n_waters, 0.0);

    ctx.list_sh.assign(ctx.n_max_inshell, std::vector<int>(ctx.n_shells, 0));
    ctx.nsort.assign(ctx.n_max_inshell, std::vector<int>(ctx.n_shells, 0));

    if (ctx.command_info.requested_gpu) {
        ctx.wshells->upload();
    }
}

void init_pshells() {
    auto& ctx = Context::instance();
    auto& atypes = ctx.atypes->cpu_data_p;
    auto& catypes = ctx.catypes->cpu_data_p;
    auto& coords_init = ctx.coords_init->cpu_data_p;
    auto* excluded = ctx.excluded->cpu_data_p;
    double mass, r2, rin2;

    ctx.heavy = std::make_unique<HostDeviceBuffer<bool>>(ctx.n_atoms, true, ctx.command_info.requested_gpu);
    auto* heavy = ctx.heavy->cpu_data_p;
    ctx.shell = std::make_unique<HostDeviceBuffer<bool>>(ctx.n_atoms, true, ctx.command_info.requested_gpu);
    auto* shell = ctx.shell->cpu_data_p;
    for (int i = 0; i < ctx.n_atoms; ++i) {
        heavy[i] = false;
        shell[i] = false;
    }
    rin2 = pow(shell_default * ctx.topo.exclusion_radius, 2);

    int n_heavy = 0, n_inshell = 0;

    for (int i = 0; i < ctx.n_atoms; i++) {
        mass = catypes[atypes[i].code - 1].m;
        if (mass < 4.0) {
            heavy[i] = false;
        } else {
            heavy[i] = true;
            n_heavy++;
        }

        if (heavy[i] && !excluded[i] && i < ctx.n_atoms_solute) {
            r2 = pow(coords_init[i].x - ctx.topo.solute_center.x, 2) + pow(coords_init[i].y - ctx.topo.solute_center.y, 2) + pow(coords_init[i].z - ctx.topo.solute_center.z, 2);
            if (r2 > rin2) {
                shell[i] = true;
                n_inshell++;
            } else {
                shell[i] = false;
            }
        }
    }

    if (ctx.command_info.requested_gpu) {
        ctx.heavy->upload();
        ctx.shell->upload();
    }

    printf("n_heavy = %d, n_inshell = %d\n", n_heavy, n_inshell);
}

// Marks heavy atoms for the shell/excluded arrays.
// Shared between switch-atom and centroid shell init paths.
static int mark_heavy_atoms(Context& ctx) {
    auto& atypes = ctx.atypes->cpu_data_p;
    auto& catypes = ctx.catypes->cpu_data_p;
    auto* heavy = ctx.heavy->cpu_data_p;
    int n_heavy = 0;
    for (int i = 0; i < ctx.n_atoms; i++) {
        double mass = catypes[atypes[i].code - 1].m;
        if (mass < 4.0) {
            heavy[i] = false;
        } else {
            heavy[i] = true;
            n_heavy++;
        }
    }
    return n_heavy;
}

static void allocate_pshell_buffers(Context& ctx) {
    ctx.heavy = std::make_unique<HostDeviceBuffer<bool>>(ctx.n_atoms, true, ctx.command_info.requested_gpu);
    auto* heavy = ctx.heavy->cpu_data_p;
    ctx.shell = std::make_unique<HostDeviceBuffer<bool>>(ctx.n_atoms, true, ctx.command_info.requested_gpu);
    auto* shell = ctx.shell->cpu_data_p;
    for (int i = 0; i < ctx.n_atoms; ++i) {
        heavy[i] = false;
        shell[i] = false;
    }
}

void init_pshells_from_charge_groups() {
    auto& ctx = Context::instance();
    auto& coords_init = ctx.coords_init->cpu_data_p;
    auto* excluded = ctx.excluded->cpu_data_p;
    double r2, rin2;
    auto& charge_groups = ctx.charge_group_config;
    const bool use_switch_atom = charge_groups.iuse_switch_atom == 1;

    allocate_pshell_buffers(ctx);
    auto* heavy = ctx.heavy->cpu_data_p;
    auto* shell = ctx.shell->cpu_data_p;
    rin2 = pow(ctx.md.shell_radius, 2);

    int n_heavy = mark_heavy_atoms(ctx);
    int n_inshell = 0;

    for (int grp = 0; grp < charge_groups.n_cgrps_solute; grp++) {
        const auto& charge_group = charge_groups.charge_groups[grp];
        int i = charge_group.iswitch - 1;
        if (heavy[i] && !excluded[i] && i < ctx.n_atoms_solute) {
            double cx = coords_init[i].x;
            double cy = coords_init[i].y;
            double cz = coords_init[i].z;
            if (!use_switch_atom) {
                cx = 0.0;
                cy = 0.0;
                cz = 0.0;
                for (int atom : charge_group.atoms) {
                    int ai = atom - 1;
                    cx += coords_init[ai].x;
                    cy += coords_init[ai].y;
                    cz += coords_init[ai].z;
                }
                double inv_atoms = 1.0 / static_cast<double>(charge_group.atoms.size());
                cx *= inv_atoms;
                cy *= inv_atoms;
                cz *= inv_atoms;
            }

            r2 = pow(cx - ctx.topo.solute_center.x, 2) + pow(cy - ctx.topo.solute_center.y, 2) + pow(cz - ctx.topo.solute_center.z, 2);
            bool in_shell = r2 > rin2;
            for (int atom : charge_group.atoms) {
                shell[atom - 1] = in_shell;
                if (in_shell) {
                    n_inshell++;
                }
            }
        }
    }

    if (ctx.command_info.requested_gpu) {
        ctx.heavy->upload();
        ctx.shell->upload();
    }

    printf("(%s): n_heavy = %d, n_inshell = %d\n", use_switch_atom ? "switch atoms" : "centroids", n_heavy, n_inshell);
}


/* =============================================
 * == CALCUTED IN THE INTEGRATION
 * =============================================
 */

void init_patoms() {
    auto& ctx = Context::instance();
    ctx.n_patoms = ctx.n_atoms_solute - ctx.n_qatoms;

    ctx.p_atoms.resize(ctx.n_patoms);

    // Loop through all solutes, adding a p atom to the list every time a non-q atom is encountered
    int pi = 0;
    int qi = 0;
    for (int i = 0; i < ctx.n_atoms_solute; i++) {
        if (qi < ctx.n_qatoms && i == ctx.q_atoms[qi]) {
            qi++;
        } else {
            ctx.p_atoms[pi] = i;
            pi++;
        }
    }
}

// Build a compact atom/state lookup layer on top of the normal and Q parameter
// tables. Base codes 1..n_atoms hold the default atom entries, while only Q
// states above 0 get appended as extra codes.
void init_unified_atom_parameters() {
    auto& ctx = Context::instance();
    auto& charges = ctx.charges->cpu_data_p;
    auto& ccharges = ctx.ccharges->cpu_data_p;
    auto& atypes = ctx.atypes->cpu_data_p;
    auto& catypes = ctx.catypes->cpu_data_p;
    const int n_states = ctx.n_parameter_states();

    ctx.atom_to_qi.assign(ctx.n_atoms, -1);
    for (int qi = 0; qi < ctx.n_qatoms; qi++) {
        ctx.atom_to_qi[ctx.q_atoms[qi]] = qi;
    }

    const int n_extra_q_states = n_states > 0 ? n_states - 1 : 0;
    const int n_codes = ctx.n_atoms + ctx.n_qatoms * n_extra_q_states;
    ctx.unified_ccharges = std::make_unique<HostDeviceBuffer<ccharge_t>>(n_codes, true, ctx.command_info.requested_gpu);
    ctx.unified_catypes = std::make_unique<HostDeviceBuffer<catype_t>>(n_codes, true, ctx.command_info.requested_gpu);
    auto* unified_ccharges = ctx.unified_ccharges->cpu_data_p;
    auto* unified_catypes = ctx.unified_catypes->cpu_data_p;

    for (int atom_idx = 0; atom_idx < ctx.n_atoms; atom_idx++) {
        const int unified_code = atom_idx + 1;
        const int qi = ctx.atom_to_qi[atom_idx];

        ccharge_t resolved_charge = {};
        resolved_charge.code = unified_code;

        catype_t resolved_type = {};
        resolved_type.code = unified_code;

        if (qi != -1 && ctx.n_lambdas > 0) {
            resolved_charge.charge = ctx.q_charges[qi].charge;
            resolved_type = ctx.q_catypes[ctx.q_atypes[qi].code - 1];
            resolved_type.code = unified_code;
        } else {
            resolved_charge.charge = ccharges[charges[atom_idx].code - 1].charge;
            resolved_type = catypes[atypes[atom_idx].code - 1];
            resolved_type.code = unified_code;
        }

        unified_ccharges[unified_code - 1] = resolved_charge;
        unified_catypes[unified_code - 1] = resolved_type;
    }

    for (int state = 1; state < n_states; state++) {
        for (int qi = 0; qi < ctx.n_qatoms; qi++) {
            const int unified_code = ctx.n_atoms + (state - 1) * ctx.n_qatoms + qi + 1;

            ccharge_t resolved_charge = {};
            resolved_charge.code = unified_code;
            resolved_charge.charge = ctx.q_charges[qi + state * ctx.n_qatoms].charge;

            catype_t resolved_type = ctx.q_catypes[ctx.q_atypes[qi + state * ctx.n_qatoms].code - 1];
            resolved_type.code = unified_code;

            unified_ccharges[unified_code - 1] = resolved_charge;
            unified_catypes[unified_code - 1] = resolved_type;
        }
    }

    if (ctx.command_info.requested_gpu) {
        ctx.unified_ccharges->upload();
        ctx.unified_catypes->upload();
    }
}
void stop_cm_translation(Context& ctx) {
    auto& atypes = ctx.atypes->cpu_data_p;
    auto& catypes = ctx.catypes->cpu_data_p;
    auto& velocities = ctx.velocities->cpu_data_p;
    double total_mass = 0;
    coord_t vcm = {};

    for (int ai = 0; ai < ctx.n_atoms; ai++) {
        const double rmass = catypes[atypes[ai].code - 1].m;
        total_mass += rmass;
        vcm.x += velocities[ai].x * rmass;
        vcm.y += velocities[ai].y * rmass;
        vcm.z += velocities[ai].z * rmass;
    }

    vcm.x /= total_mass;
    vcm.y /= total_mass;
    vcm.z /= total_mass;

    for (int ai = 0; ai < ctx.n_atoms; ai++) {
        velocities[ai].x -= vcm.x;
        velocities[ai].y -= vcm.y;
        velocities[ai].z -= vcm.z;
    }
}
