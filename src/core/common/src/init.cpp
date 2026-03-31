#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>
#include <unistd.h>

#include "constants.h"
#include "context.h"
#include "cpu_handler.h"
#include "cpu_shake.h"
#include "cpu_utils.h"
#include "output.h"
#include "parse.h"
#include "init.h"

// Remove bonds, angles, torsions and impropers which are excluded or changed in the FEP file
void exclude_qatom_definitions() {
    auto& ctx = Context::instance();
    int excluded;
    int ai = 0, bi = 0, ii = 0, ti = 0;
    int qai = 0, qbi = 0, qii = 0, qti = 0;

    excluded = 0;
    int solute_excluded = 0;
    if (ctx.n_qangles > 0) {
        for (int i = 0; i < ctx.n_angles; i++) {
            if (ctx.angles[i].ai == ctx.q_angles[qai].ai && ctx.angles[i].aj == ctx.q_angles[qai].aj && ctx.angles[i].ak == ctx.q_angles[qai].ak) {
                qai++;
                excluded++;
                if (i < ctx.n_angles_solute) solute_excluded++;
            } else {
                ctx.angles[ai] = ctx.angles[i];
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
            if (ctx.bonds[i].ai == ctx.q_bonds[qbi].ai && ctx.bonds[i].aj == ctx.q_bonds[qbi].aj) {
                qbi++;
                excluded++;
                if (i < ctx.n_bonds_solute) solute_excluded++;
            } else {
                ctx.bonds[bi] = ctx.bonds[i];
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
            if (ctx.torsions[i].ai == ctx.q_torsions[qti].ai && ctx.torsions[i].aj == ctx.q_torsions[qti].aj && ctx.torsions[i].ak == ctx.q_torsions[qti].ak && ctx.torsions[i].al == ctx.q_torsions[qti].al) {
                qti++;
                excluded++;
                if (i < ctx.n_torsions_solute) solute_excluded++;
            } else {
                ctx.torsions[ti] = ctx.torsions[i];
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
    int n_excl;
    int ai = 0, bi = 0, ii = 0, ti = 0;

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
        if (ctx.excluded[ctx.impropers[i].ai - 1] && ctx.excluded[ctx.impropers[i].aj - 1] && ctx.excluded[ctx.impropers[i].ak - 1] && ctx.excluded[ctx.impropers[i].al - 1]) {
            n_excl++;
        } else {
            ctx.impropers[ii] = ctx.impropers[i];
            ii++;
        }
    }
    printf("original: %d. # excluded impropers: %d\n", ctx.n_impropers, n_excl);
    ctx.n_impropers -= n_excl;

    n_excl = 0;
    for (int i = 0; i < ctx.n_torsions; i++) {
        if (ctx.excluded[ctx.torsions[i].ai - 1] && ctx.excluded[ctx.torsions[i].aj - 1] && ctx.excluded[ctx.torsions[i].ak - 1] && ctx.excluded[ctx.torsions[i].al - 1]) {
            n_excl++;
        } else {
            ctx.torsions[ti] = ctx.torsions[i];
            ti++;
        }
    }
    printf("original: %d. # excluded torsions: %d\n", ctx.n_torsions, n_excl);
    ctx.n_torsions -= n_excl;
}

void exclude_shaken_definitions() {
    auto& ctx = Context::instance();
    int excluded;
    int solute_excluded;
    int bi = 0;
    int si = 0;
    int ang_i = 0;
    int ai, aj;

    excluded = 0;
    solute_excluded = 0;
    if (ctx.n_shake_constraints > 0) {
        for (int i = 0; i < ctx.n_bonds; i++) {
            if (ctx.bonds[i].ai == ctx.shake_bonds[si].ai && ctx.bonds[i].aj == ctx.shake_bonds[si].aj) {
                si++;
                excluded++;
                if (i < ctx.n_bonds_solute) solute_excluded++;
            } else {
                ctx.bonds[bi] = ctx.bonds[i];
                bi++;
            }
        }
        ctx.n_bonds -= excluded;
        ctx.n_bonds_solute -= solute_excluded;
    }

    excluded = 0;
    solute_excluded = 0;
    if (ctx.n_shake_constraints > 0) {
        // Mark angles whose terminal atoms (i and k) match a shaken bond
        for (int i = 0; i < ctx.n_shake_constraints; i++) {
            ai = ctx.shake_bonds[i].ai;
            aj = ctx.shake_bonds[i].aj;
            for (int j = 0; j < ctx.n_angles; j++) {
                if ((ctx.angles[j].ai == ai && ctx.angles[j].ak == aj) || (ctx.angles[j].ai == aj && ctx.angles[j].ak == ai)) {
                    ctx.angles[j].code = 0;
                    break;
                }
            }
        }

        for (int i = 0; i < ctx.n_angles; i++) {
            if (ctx.angles[i].code == 0) {
                excluded++;
                if (i < ctx.n_angles_solute) solute_excluded++;
            } else {
                ctx.angles[ang_i] = ctx.angles[i];
                ang_i++;
            }
        }
        ctx.n_angles -= excluded;
        ctx.n_angles_solute -= solute_excluded;
    }
}

void init_velocities() {
    auto& ctx = Context::instance();
    ctx.velocities.resize(ctx.n_atoms);

    // If not previous value set, use a Maxwell distribution to fill velocities
    double kT = Boltz * ctx.md.initial_temperature;
    double sd, mass;
    for (int i = 0; i < ctx.n_atoms; i++) {
        mass = ctx.catypes[ctx.atypes[i].code - 1].m;
        sd = sqrt(kT / mass);

        ctx.velocities[i].x = gauss(0, sd);
        ctx.velocities[i].y = gauss(0, sd);
        ctx.velocities[i].z = gauss(0, sd);
    }
}

void init_dvelocities() {
    auto& ctx = Context::instance();
    ctx.dvelocities.assign(ctx.n_atoms, {});
}

void init_xcoords() {
    auto& ctx = Context::instance();
    ctx.xcoords.resize(ctx.n_atoms);
}

void init_inv_mass() {
    auto& ctx = Context::instance();
    ctx.winv.resize(ctx.n_atoms);
    for (int ai = 0; ai < ctx.n_atoms; ai++) {
        ctx.winv[ai] = 1 / ctx.catypes[ctx.atypes[ai].code - 1].m;
    }
}

/* =============================================
 * == BOUNDARY RESTRAINTS
 * =============================================
 */

void init_water_sphere() {
    auto& ctx = Context::instance();
    ctx.Dwmz = 0.26 * exp(-0.19 * (ctx.topo.solvent_radius - 15)) + 0.74;
    ctx.awmz = 0.2 / (1 + exp(0.4 * (ctx.topo.solvent_radius - 25))) + 0.3;

    printf("Dwmz = %f, awmz = %f\n", ctx.Dwmz, ctx.awmz);
}

// ONLY call if there are actually solvent atoms, or get segfaulted
void init_wshells() {
    auto& ctx = Context::instance();
    int n_inshell;
    double drs, router, ri, dr, Vshell, rshell;
    if (ctx.mu_w == 0) {
        // Get water properties from first water molecule
        cbond_t cbondw = ctx.cbonds[ctx.bonds[ctx.n_atoms_solute].code - 1];
        cangle_t canglew = ctx.cangles[ctx.angles[ctx.n_atoms_solute].code - 1];

        ccharge_t ccharge_ow = ctx.ccharges[ctx.charges[ctx.n_atoms_solute].code - 1];
        ctx.crg_ow = ccharge_ow.charge;

        ctx.mu_w = -ctx.crg_ow * cbondw.b0 * cos(canglew.th0 / 2);
    }

    drs = wpolr_layer / drouter;

    ctx.n_shells = (int)floor(-0.5 + sqrt(2 * drs + 0.25));
    ctx.wshells.resize(ctx.n_shells);

    printf("n_shells = %d\n", ctx.n_shells);

    router = ctx.topo.solvent_radius;
    ctx.n_max_inshell = 0;

    for (int i = 0; i < ctx.n_shells; i++) {
        ctx.wshells[i].avtheta = 0;
        ctx.wshells[i].avn_inshell = 0;
        ctx.wshells[i].router = router;
        dr = drouter * (i + 1);
        ri = router - dr;
        ctx.wshells[i].dr = dr;
        Vshell = pow(router, 3) - pow(ri, 3);
        n_inshell = (int)floor(4 * M_PI / 3 * Vshell * rho_water);
        if (n_inshell > ctx.n_max_inshell) {
            ctx.n_max_inshell = n_inshell;
        }
        rshell = pow(0.5 * (pow(router, 3) + pow(ri, 3)), 1.0 / 3.0);

        // --- Note below: 0.98750 = (1-1/epsilon) for water
        ctx.wshells[i].cstb = ctx.crgQtot * 0.98750 / (rho_water * ctx.mu_w * 4 * M_PI * pow(rshell, 2));

        router -= dr;
    }

    // rc > wshells[n_shells-1].router - wshells[n_shells-1].dr
    printf("shell 0: (%f, %f). shell 1: (%f, %f). shell 2: (%f, %f).\n", ctx.wshells[0].router, ctx.wshells[0].router - ctx.wshells[0].dr, ctx.wshells[1].router, ctx.wshells[1].router - ctx.wshells[1].dr, ctx.wshells[2].router, ctx.wshells[2].router - ctx.wshells[2].dr);

    ctx.n_max_inshell = ctx.n_waters;  // Make largest a little bigger just in case

    // Initialize arrays needed for bookkeeping
    ctx.theta.assign(ctx.n_waters, 0.0);
    ctx.theta0.assign(ctx.n_waters, 0.0);
    ctx.tdum.assign(ctx.n_waters, 0.0);

    ctx.list_sh.assign(ctx.n_max_inshell, std::vector<int>(ctx.n_shells, 0));
    ctx.nsort.assign(ctx.n_max_inshell, std::vector<int>(ctx.n_shells, 0));
}

void init_pshells() {
    auto& ctx = Context::instance();
    double mass, r2, rin2;

    ctx.heavy = std::make_unique<bool[]>(ctx.n_atoms);
    ctx.shell = std::make_unique<bool[]>(ctx.n_atoms);
    for (int i = 0; i < ctx.n_atoms; ++i) {
        ctx.heavy[i] = false;
        ctx.shell[i] = false;
    }
    rin2 = pow(shell_default * ctx.topo.exclusion_radius, 2);

    int n_heavy = 0, n_inshell = 0;

    for (int i = 0; i < ctx.n_atoms; i++) {
        mass = ctx.catypes[ctx.atypes[i].code - 1].m;
        if (mass < 4.0) {
            ctx.heavy[i] = false;
        } else {
            ctx.heavy[i] = true;
            n_heavy++;
        }

        if (ctx.heavy[i] && !ctx.excluded[i] && i < ctx.n_atoms_solute) {
            r2 = pow(ctx.coords_top[i].x - ctx.topo.solute_center.x, 2) + pow(ctx.coords_top[i].y - ctx.topo.solute_center.y, 2) + pow(ctx.coords_top[i].z - ctx.topo.solute_center.z, 2);
            if (r2 > rin2) {
                ctx.shell[i] = true;
                n_inshell++;
            } else {
                ctx.shell[i] = false;
            }
        }
    }

    printf("n_heavy = %d, n_inshell = %d\n", n_heavy, n_inshell);
}

// Marks heavy atoms for the shell/excluded arrays.
// Shared between switch-atom and centroid shell init paths.
static int mark_heavy_atoms(Context& ctx) {
    int n_heavy = 0;
    for (int i = 0; i < ctx.n_atoms; i++) {
        double mass = ctx.catypes[ctx.atypes[i].code - 1].m;
        if (mass < 4.0) {
            ctx.heavy[i] = false;
        } else {
            ctx.heavy[i] = true;
            n_heavy++;
        }
    }
    return n_heavy;
}

// Shell init using switch atom distance (matches Fortran make_shell).
// Used when iuse_switch_atom == 1.
void init_pshells_with_switch_atoms() {
    auto& ctx = Context::instance();
    double r2, rin2;

    ctx.heavy = std::make_unique<bool[]>(ctx.n_atoms);
    ctx.shell = std::make_unique<bool[]>(ctx.n_atoms);
    for (int i = 0; i < ctx.n_atoms; ++i) {
        ctx.heavy[i] = false;
        ctx.shell[i] = false;
    }
    rin2 = pow(ctx.md.shell_radius, 2);

    int n_heavy = mark_heavy_atoms(ctx);
    int n_inshell = 0;

    for (int grp = 0; grp < ctx.n_cgrps_solute; grp++) {
        cgrp_t cgrp = ctx.charge_groups[grp];
        int i = cgrp.iswitch - 1;
        if (ctx.heavy[i] && !ctx.excluded[i] && i < ctx.n_atoms_solute) {
            r2 = pow(ctx.coords_top[i].x - ctx.topo.solute_center.x, 2) + pow(ctx.coords_top[i].y - ctx.topo.solute_center.y, 2) + pow(ctx.coords_top[i].z - ctx.topo.solute_center.z, 2);
            bool in_shell = r2 > rin2;
            for (int j = 0; j < cgrp.n_atoms; j++) {
                ctx.shell[cgrp.a[j] - 1] = in_shell;
                if (in_shell) {
                    n_inshell++;
                }
            }
        }
    }

    printf("(switch atoms): n_heavy = %d, n_inshell = %d\n", n_heavy, n_inshell);
}

// Shell init using charge group centroid distance (matches Fortran make_shell2).
// Used when iuse_switch_atom == 0.
void init_pshells_with_centroids() {
    auto& ctx = Context::instance();
    double r2, rin2;

    ctx.heavy = std::make_unique<bool[]>(ctx.n_atoms);
    ctx.shell = std::make_unique<bool[]>(ctx.n_atoms);
    for (int i = 0; i < ctx.n_atoms; ++i) {
        ctx.heavy[i] = false;
        ctx.shell[i] = false;
    }
    rin2 = pow(ctx.md.shell_radius, 2);

    int n_heavy = mark_heavy_atoms(ctx);
    int n_inshell = 0;

    for (int grp = 0; grp < ctx.n_cgrps_solute; grp++) {
        cgrp_t cgrp = ctx.charge_groups[grp];
        int i = cgrp.iswitch - 1;
        if (ctx.heavy[i] && !ctx.excluded[i] && i < ctx.n_atoms_solute) {
            // Compute centroid of charge group
            double cx = 0, cy = 0, cz = 0;
            for (int j = 0; j < cgrp.n_atoms; j++) {
                int ai = cgrp.a[j] - 1;
                cx += ctx.coords_top[ai].x;
                cy += ctx.coords_top[ai].y;
                cz += ctx.coords_top[ai].z;
            }
            cx /= cgrp.n_atoms;
            cy /= cgrp.n_atoms;
            cz /= cgrp.n_atoms;

            r2 = pow(cx - ctx.topo.solute_center.x, 2) + pow(cy - ctx.topo.solute_center.y, 2) + pow(cz - ctx.topo.solute_center.z, 2);
            bool in_shell = r2 > rin2;
            for (int j = 0; j < cgrp.n_atoms; j++) {
                ctx.shell[cgrp.a[j] - 1] = in_shell;
                if (in_shell) {
                    n_inshell++;
                }
            }
        }
    }

    printf("(centroids): n_heavy = %d, n_inshell = %d\n", n_heavy, n_inshell);
}

void init_restrseqs() {
    auto& ctx = Context::instance();
    ctx.n_restrseqs = 1;
    ctx.restrseqs.resize(1);

    restrseq_t seq;
    seq.ai = 1;
    seq.aj = 14;
    seq.k = 1.0;
    seq.ih = 0;
    seq.to_center = 2;

    ctx.restrseqs[0] = seq;
}

/* =============================================
 * == SHAKE
 * =============================================
 */

void init_shake() {
    auto& ctx = Context::instance();
    int ai, aj;
    int mol = 0;
    int shake;
    int n_solute_shake_constraints = 0;
    double excl_shake = 0;

    ctx.n_shake_constraints = 0;
    ctx.mol_n_shakes.assign(ctx.n_molecules, 0);

    for (int bi = 0; bi < ctx.n_bonds; bi++) {
        ai = ctx.bonds[bi].ai - 1;
        aj = ctx.bonds[bi].aj - 1;

        while (mol + 1 < ctx.n_molecules && ai + 1 >= ctx.molecules[mol + 1]) {
            // new molecule
            mol += 1;
        }

        if ((ctx.md.shake_hydrogens && (!ctx.heavy[ai] || !ctx.heavy[aj])) || (ctx.md.shake_solute && ai + 1 <= ctx.n_atoms_solute) || (ctx.md.shake_solvent && ai + 1 > ctx.n_atoms_solute)) {
            ctx.mol_n_shakes[mol]++;
            ctx.n_shake_constraints++;

            if (ctx.excluded[ai]) excl_shake += 0.5;
            if (ctx.excluded[aj]) excl_shake += 0.5;
        }
    }

    ctx.shake_bonds.resize(ctx.n_shake_constraints);
    mol = 0;
    shake = 0;
    for (int bi = 0; bi < ctx.n_bonds; bi++) {
        ai = ctx.bonds[bi].ai - 1;
        aj = ctx.bonds[bi].aj - 1;

        while (mol + 1 < ctx.n_molecules && ai + 1 >= ctx.molecules[mol + 1]) {
            // new molecule
            mol += 1;
        }

        if ((ctx.md.shake_hydrogens && (!ctx.heavy[ai] || !ctx.heavy[aj])) || (ctx.md.shake_solute && ai + 1 <= ctx.n_atoms_solute) || (ctx.md.shake_solvent && ai + 1 > ctx.n_atoms_solute)) {
            ctx.shake_bonds[shake].ai = ai + 1;
            ctx.shake_bonds[shake].aj = aj + 1;
            ctx.shake_bonds[shake].dist2 = pow(ctx.cbonds[ctx.bonds[bi].code - 1].b0, 2);
            shake++;
        }
    }

    // Get total number of shake constraints in solute (used for separate scaling of temperatures)
    for (int i = 0; i < ctx.n_molecules - ctx.n_waters; i++) {
        n_solute_shake_constraints += ctx.mol_n_shakes[i];
    }

    ctx.Ndegf = 3 * ctx.n_atoms - ctx.n_shake_constraints;
    ctx.Ndegfree = ctx.Ndegf - 3 * ctx.n_excluded + excl_shake;

    ctx.Ndegf_solvent = ctx.Ndegf - 3 * ctx.n_atoms_solute + n_solute_shake_constraints;
    ctx.Ndegf_solute = ctx.Ndegf - ctx.Ndegf_solvent;

    ctx.Ndegfree_solvent = ctx.Ndegfree - (ctx.n_shake_constraints - n_solute_shake_constraints);
    ctx.Ndegfree_solute = ctx.Ndegfree - ctx.Ndegfree_solvent;

    printf("n_shake_constrains = %d, n_solute_shake_constraints = %d, excl_shake = %f\n", ctx.n_shake_constraints, n_solute_shake_constraints, excl_shake);

    if (ctx.Ndegfree_solvent * ctx.Ndegfree_solute == 0) {
        ctx.separate_scaling = false;
    } else {
        ctx.separate_scaling = true;
    }
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
        if (ctx.n_qatoms > 0 && qi < ctx.n_qatoms && i == ctx.q_atoms[qi].a - 1) {
            qi++;
        } else {
            ctx.p_atoms[pi].a = i + 1;
            pi++;
        }
    }
}

void init_variables() {
    auto& ctx = Context::instance();
    // From MD file
    init_md("md.csv");

    ctx.dt = time_unit * ctx.md.stepsize;
    ctx.tau_T = time_unit * ctx.md.bath_coupling;

    if (ctx.run_gpu && ctx.n_lambdas > 2) {
        printf(">>> FATAL: More than 2 states not supported on GPU architecture. Exiting...\n");
        exit(EXIT_FAILURE);
    }

    // From topology file
    init_topo("topo.csv");

    init_angles("angles.csv");
    init_atypes("atypes.csv");
    init_bonds("bonds.csv");
    init_cangles("cangles.csv");
    init_catypes("catypes.csv");
    init_cbonds("cbonds.csv");
    init_ccharges("ccharges.csv");
    init_charges("charges.csv");
    init_cimpropers("cimpropers.csv");
    init_coords("coords.csv");
    init_ctorsions("ctorsions.csv");
    init_excluded("excluded.csv");
    init_molecules("molecules.csv");
    init_impropers("impropers.csv");
    init_torsions("torsions.csv");
    init_LJ_matrix();
    init_ngbrs14("ngbrs14.csv");
    init_ngbrs23("ngbrs23.csv");
    init_ngbrs14_long("ngbrs14long.csv");
    init_ngbrs23_long("ngbrs23long.csv");
    // init_restrseqs();
    init_inv_mass();
    if (ctx.md.charge_groups) {
        init_charge_groups("charge_groups.csv");
    }

    // From FEP file
    init_qangcouples("q_angcouples.csv");
    init_qatoms("q_atoms.csv");
    init_qcangles("q_cangles.csv");
    init_qcatypes("q_catypes.csv");
    init_qcbonds("q_cbonds.csv");
    init_qcimpropers("q_cimpropers.csv");
    init_qctorsions("q_ctorsions.csv");
    init_qoffdiags("q_offdiags.csv");
    init_qimprcouples("q_imprcouples.csv");
    init_qsoftpairs("q_softpairs.csv");
    init_qtorcouples("q_torcouples.csv");

    init_qangles("q_angles.csv");
    init_qatypes("q_atypes.csv");
    init_qbonds("q_bonds.csv");
    init_qcharges("q_charges.csv");
    init_qelscales("q_elscales.csv");
    init_qexclpairs("q_exclpairs.csv");
    init_qimpropers("q_impropers.csv");
    init_qshakes("q_shakes.csv");
    init_qsoftcores("q_softcores.csv");
    init_qtorsions("q_torsions.csv");

    // First part of shrink topology, this needs to be done first as shake constraints are based on bonds
    exclude_qatom_definitions();
    // exclude_all_atoms_excluded_definitions();

    // Shake constraints, need to be initialized before last part of shrink_topology
    if (ctx.md.charge_groups) {
        if (ctx.iuse_switch_atom == 1) {
            init_pshells_with_switch_atoms();
        } else {
            init_pshells_with_centroids();
        }
    } else {
        init_pshells();
    }
    init_shake();

    // Now remove shaken bonds
    exclude_shaken_definitions();

    // Init random seed from MD file
    srand(ctx.md.random_seed);

    // From calculation in the integration
    init_patoms();
    init_velocities();
    init_dvelocities();
    init_xcoords();

    // From input file
    init_icoords("i_coords.csv");
    init_ivelocities("i_velocities.csv");

    // Init waters, boundary restrains
    ctx.n_waters = (ctx.n_atoms - ctx.n_atoms_solute) / 3;
    if (ctx.n_waters > 0) {
        init_water_sphere();
        init_wshells();
    }

    // Init energy
    ctx.EQ_total.assign(ctx.n_lambdas, {});
    ctx.EQ_bond.assign(ctx.n_lambdas, {});
    ctx.EQ_nonbond_qq.assign(ctx.n_lambdas, {});
    ctx.EQ_nonbond_qp.assign(ctx.n_lambdas, {});
    ctx.EQ_nonbond_qw.assign(ctx.n_lambdas, {});
    ctx.EQ_nonbond_qx.assign(ctx.n_lambdas, {});
    ctx.EQ_restraint.assign(ctx.n_lambdas, {});

    if (ctx.n_shake_constraints > 0) {
        initial_shaking();
        stop_cm_translation();
    }

    // Write header to file
    write_header("coords.csv");
    write_header("velocities.csv");
    write_energy_header();
}

void clean_variables() {
    auto& ctx = Context::instance();

    for (auto& charge_group : ctx.charge_groups) {
        delete[] charge_group.a;
        charge_group.a = nullptr;
    }
    ctx.charge_groups.clear();

    ctx.angles.clear();
    ctx.atypes.clear();
    ctx.bonds.clear();
    ctx.cangles.clear();
    ctx.catypes.clear();
    ctx.cbonds.clear();
    ctx.ccharges.clear();
    ctx.charges.clear();
    ctx.cimpropers.clear();
    ctx.coords.clear();
    ctx.coords_top.clear();
    ctx.ctorsions.clear();
    ctx.excluded.reset();
    ctx.heavy.reset();
    ctx.impropers.clear();
    ctx.torsions.clear();
    ctx.LJ_matrix.clear();
    ctx.molecules.clear();
    ctx.q_angcouples.clear();
    ctx.q_atoms.clear();
    ctx.q_cangles.clear();
    ctx.q_catypes.clear();
    ctx.q_cbonds.clear();
    ctx.q_cimpropers.clear();
    ctx.q_ctorsions.clear();
    ctx.q_offdiags.clear();
    ctx.q_imprcouples.clear();
    ctx.q_softpairs.clear();
    ctx.q_torcouples.clear();
    ctx.q_atypes.clear();
    ctx.q_charges.clear();
    ctx.q_angles.clear();
    ctx.q_bonds.clear();
    ctx.q_elscales.clear();
    ctx.q_exclpairs.clear();
    ctx.q_impropers.clear();
    ctx.q_shakes.clear();
    ctx.q_softcores.clear();
    ctx.q_torsions.clear();
    ctx.wshells.clear();
    ctx.theta.clear();
    ctx.theta0.clear();
    ctx.tdum.clear();
    ctx.list_sh.clear();
    ctx.nsort.clear();
    ctx.restrseqs.clear();
    ctx.shell.reset();
    ctx.velocities.clear();
    ctx.dvelocities.clear();
    ctx.xcoords.clear();
    ctx.mol_n_shakes.clear();
    ctx.shake_bonds.clear();
    ctx.EQ_total.clear();
    ctx.EQ_bond.clear();
    ctx.EQ_nonbond_qq.clear();
    ctx.EQ_nonbond_qp.clear();
    ctx.EQ_nonbond_qw.clear();
    ctx.EQ_nonbond_qx.clear();
    ctx.EQ_restraint.clear();
}
