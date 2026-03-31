#include "cpu/include/context.h"
#include "cpu/include/constants.h"
#include "system.h"
#include "utils.h"
#include "parse.h"
#include "bonded.h"
#include "nonbonded.h"
#include "solvent.h"
#include "restraints.h"
#include "qatoms.h"
#include "shake.h"
#include "cuda/include/CudaHandler.cuh"


#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <time.h>
#include <unistd.h>


// Remove bonds, angles, torsions and impropers which are excluded or changed in the FEP file
void exclude_qatom_definitions() {
    auto &ctx = Context::instance();
    int excluded;
    int ai = 0, bi = 0, ii = 0, ti = 0;
    int qai = 0, qbi = 0, qii = 0, qti = 0;

    excluded = 0;
    int solute_excluded = 0;
    if (ctx.n_qangles > 0) {
        for (int i = 0; i < ctx.n_angles; i++) {
            if (ctx.angles[i].ai == ctx.q_angles[qai].ai
             && ctx.angles[i].aj == ctx.q_angles[qai].aj
             && ctx.angles[i].ak == ctx.q_angles[qai].ak) {
                qai++;
                excluded++;
                if (i < ctx.n_angles_solute) solute_excluded++;
            }
            else {
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
            if (ctx.bonds[i].ai == ctx.q_bonds[qbi].ai
             && ctx.bonds[i].aj == ctx.q_bonds[qbi].aj) {
                qbi++;
                excluded++;
                if (i < ctx.n_bonds_solute) solute_excluded++;
            }
            else {
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
            if (ctx.torsions[i].ai == ctx.q_torsions[qti].ai
             && ctx.torsions[i].aj == ctx.q_torsions[qti].aj
             && ctx.torsions[i].ak == ctx.q_torsions[qti].ak
             && ctx.torsions[i].al == ctx.q_torsions[qti].al) {
                qti++;
                excluded++;
                if (i < ctx.n_torsions_solute) solute_excluded++;
            }
            else {
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
    auto &ctx = Context::instance();
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
        if (ctx.excluded[ctx.impropers[i].ai - 1]
            && ctx.excluded[ctx.impropers[i].aj - 1]
            && ctx.excluded[ctx.impropers[i].ak - 1]
            && ctx.excluded[ctx.impropers[i].al - 1]) {
            n_excl++;
        }
        else {
            ctx.impropers[ii] = ctx.impropers[i];
            ii++;
        }
    }
    printf("original: %d. # excluded impropers: %d\n", ctx.n_impropers, n_excl);
    ctx.n_impropers -= n_excl;

    n_excl = 0;
    for (int i = 0; i < ctx.n_torsions; i++) {
        if (ctx.excluded[ctx.torsions[i].ai - 1]
            && ctx.excluded[ctx.torsions[i].aj - 1]
            && ctx.excluded[ctx.torsions[i].ak - 1]
            && ctx.excluded[ctx.torsions[i].al - 1]) {
            n_excl++;
        }
        else {
            ctx.torsions[ti] = ctx.torsions[i];
            ti++;
        }
    }
    printf("original: %d. # excluded torsions: %d\n", ctx.n_torsions, n_excl);
    ctx.n_torsions -= n_excl;
}

void exclude_shaken_definitions() {
    auto &ctx = Context::instance();
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
            if (ctx.bonds[i].ai == ctx.shake_bonds[si].ai
             && ctx.bonds[i].aj == ctx.shake_bonds[si].aj) {
                si++;
                excluded++;
                if (i < ctx.n_bonds_solute) solute_excluded++;
            }
            else {
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
                if ( (ctx.angles[j].ai == ai && ctx.angles[j].ak == aj)
                    || (ctx.angles[j].ai == aj && ctx.angles[j].ak == ai) ) {
                    ctx.angles[j].code = 0;
                    break;
                }
            }
        }

        for (int i = 0; i < ctx.n_angles; i++) {
            if (ctx.angles[i].code == 0) {
                excluded++;
                if (i < ctx.n_angles_solute) solute_excluded++;
            }
            else {
                ctx.angles[ang_i] = ctx.angles[i];
                ang_i++;
            }
        }
        ctx.n_angles -= excluded;
        ctx.n_angles_solute -= solute_excluded;
    }
}

void init_velocities() {
    auto &ctx = Context::instance();
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
    auto &ctx = Context::instance();
    ctx.dvelocities.assign(ctx.n_atoms, {});
}

void init_xcoords() {
    auto &ctx = Context::instance();
    ctx.xcoords.resize(ctx.n_atoms);
}

void init_inv_mass() {
    auto &ctx = Context::instance();
    ctx.winv.resize(ctx.n_atoms);
    for (int ai = 0; ai < ctx.n_atoms; ai++) {
        ctx.winv[ai] = 1 / ctx.catypes[ctx.atypes[ai].code-1].m;
    }
}

/* =============================================
 * == BOUNDARY RESTRAINTS
 * =============================================
 */

void init_water_sphere() {
    auto &ctx = Context::instance();
    ctx.Dwmz = 0.26 * exp(-0.19 * (ctx.topo.solvent_radius - 15)) + 0.74;
    ctx.awmz = 0.2 / (1 + exp(0.4 * (ctx.topo.solvent_radius - 25))) + 0.3;

    printf("Dwmz = %f, awmz = %f\n", ctx.Dwmz, ctx.awmz);
}

//ONLY call if there are actually solvent atoms, or get segfaulted
void init_wshells() {
    auto &ctx = Context::instance();
    int n_inshell;
    double drs, router, ri, dr, Vshell, rshell;
    if (ctx.mu_w == 0) {
        // Get water properties from first water molecule
        cbond_t cbondw = ctx.cbonds[ctx.bonds[ctx.n_atoms_solute].code-1];
        cangle_t canglew = ctx.cangles[ctx.angles[ctx.n_atoms_solute].code-1];

        ccharge_t ccharge_ow = ctx.ccharges[ctx.charges[ctx.n_atoms_solute].code - 1];
        ctx.crg_ow = ccharge_ow.charge;

        ctx.mu_w = -ctx.crg_ow * cbondw.b0 * cos(canglew.th0 / 2);
    }

    drs = wpolr_layer / drouter;

    ctx.n_shells = (int) floor(-0.5 + sqrt(2*drs + 0.25));
    ctx.wshells.resize(ctx.n_shells);

    printf("n_shells = %d\n", ctx.n_shells);

    router = ctx.topo.solvent_radius;
    ctx.n_max_inshell = 0;

    for (int i = 0; i < ctx.n_shells; i++) {
        ctx.wshells[i].avtheta = 0;
        ctx.wshells[i].avn_inshell = 0;
        ctx.wshells[i].router = router;
        dr = drouter * (i+1);
        ri = router - dr;
        ctx.wshells[i].dr = dr;
        Vshell = pow(router, 3) - pow(ri, 3);
        n_inshell = (int) floor(4 * M_PI / 3 * Vshell * rho_water);
        if (n_inshell > ctx.n_max_inshell) {
            ctx.n_max_inshell = n_inshell;
        }
        rshell = pow(0.5 * (pow(router, 3) + pow(ri, 3)), 1.0/3.0);

        // --- Note below: 0.98750 = (1-1/epsilon) for water
        ctx.wshells[i].cstb = ctx.crgQtot * 0.98750 / (rho_water * ctx.mu_w * 4 * M_PI * pow(rshell, 2));

        router -= dr;
    }

        // rc > wshells[n_shells-1].router - wshells[n_shells-1].dr
        printf("shell 0: (%f, %f). shell 1: (%f, %f). shell 2: (%f, %f).\n"
            , ctx.wshells[0].router, ctx.wshells[0].router - ctx.wshells[0].dr
            , ctx.wshells[1].router, ctx.wshells[1].router - ctx.wshells[1].dr
            , ctx.wshells[2].router, ctx.wshells[2].router - ctx.wshells[2].dr
        );

    ctx.n_max_inshell = ctx.n_waters; // Make largest a little bigger just in case

    // Initialize arrays needed for bookkeeping
    ctx.theta.assign(ctx.n_waters, 0.0);
    ctx.theta0.assign(ctx.n_waters, 0.0);
    ctx.tdum.assign(ctx.n_waters, 0.0);

    ctx.list_sh.assign(ctx.n_max_inshell, std::vector<int>(ctx.n_shells, 0));
    ctx.nsort.assign(ctx.n_max_inshell, std::vector<int>(ctx.n_shells, 0));
}

void init_pshells() {
    auto &ctx = Context::instance();
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
        mass = ctx.catypes[ctx.atypes[i].code-1].m;
        if (mass < 4.0) {
            ctx.heavy[i] = false;
        }
        else {
            ctx.heavy[i] = true;
            n_heavy++;
        }

        if (ctx.heavy[i] && !ctx.excluded[i] && i < ctx.n_atoms_solute) {
            r2 = pow(ctx.coords_top[i].x - ctx.topo.solute_center.x, 2) 
                + pow(ctx.coords_top[i].y - ctx.topo.solute_center.y, 2)
                + pow(ctx.coords_top[i].z - ctx.topo.solute_center.z, 2);
            if (r2 > rin2) {
                ctx.shell[i] = true;
                n_inshell++;
            }
            else {
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
        double mass = ctx.catypes[ctx.atypes[i].code-1].m;
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
    auto &ctx = Context::instance();
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
        int i = cgrp.iswitch-1;
        if (ctx.heavy[i] && !ctx.excluded[i] && i < ctx.n_atoms_solute) {
            r2 = pow(ctx.coords_top[i].x - ctx.topo.solute_center.x, 2)
                + pow(ctx.coords_top[i].y - ctx.topo.solute_center.y, 2)
                + pow(ctx.coords_top[i].z - ctx.topo.solute_center.z, 2);
            bool in_shell = r2 > rin2;
            for (int j = 0; j < cgrp.n_atoms; j++) {
                ctx.shell[cgrp.a[j]-1] = in_shell;
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
    auto &ctx = Context::instance();
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
        int i = cgrp.iswitch-1;
        if (ctx.heavy[i] && !ctx.excluded[i] && i < ctx.n_atoms_solute) {
            // Compute centroid of charge group
            double cx = 0, cy = 0, cz = 0;
            for (int j = 0; j < cgrp.n_atoms; j++) {
                int ai = cgrp.a[j]-1;
                cx += ctx.coords_top[ai].x;
                cy += ctx.coords_top[ai].y;
                cz += ctx.coords_top[ai].z;
            }
            cx /= cgrp.n_atoms;
            cy /= cgrp.n_atoms;
            cz /= cgrp.n_atoms;

            r2 = pow(cx - ctx.topo.solute_center.x, 2)
                + pow(cy - ctx.topo.solute_center.y, 2)
                + pow(cz - ctx.topo.solute_center.z, 2);
            bool in_shell = r2 > rin2;
            for (int j = 0; j < cgrp.n_atoms; j++) {
                ctx.shell[cgrp.a[j]-1] = in_shell;
                if (in_shell) {
                    n_inshell++;
                }
            }
        }
    }

    printf("(centroids): n_heavy = %d, n_inshell = %d\n", n_heavy, n_inshell);
}

void init_restrseqs() {
    auto &ctx = Context::instance();
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
    auto &ctx = Context::instance();
    int ai, aj;
    int mol = 0;
    int shake;
    int n_solute_shake_constraints = 0;
    double excl_shake = 0;

    ctx.n_shake_constraints = 0;
    ctx.mol_n_shakes.assign(ctx.n_molecules, 0);
    
    for (int bi = 0; bi < ctx.n_bonds; bi++) {
        ai = ctx.bonds[bi].ai-1;
        aj = ctx.bonds[bi].aj-1;

        while(mol+1 < ctx.n_molecules && ai+1 >= ctx.molecules[mol+1]) {
            // new molecule
            mol += 1;
        }

        if ( (ctx.md.shake_hydrogens && (!ctx.heavy[ai] || !ctx.heavy[aj]))
            || (ctx.md.shake_solute && ai+1 <= ctx.n_atoms_solute) 
            || (ctx.md.shake_solvent && ai+1 > ctx.n_atoms_solute) ) {
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
        ai = ctx.bonds[bi].ai-1;
        aj = ctx.bonds[bi].aj-1;

        while(mol+1 < ctx.n_molecules && ai+1 >= ctx.molecules[mol+1]) {
            // new molecule
            mol += 1;
        }

        if ( (ctx.md.shake_hydrogens && (!ctx.heavy[ai] || !ctx.heavy[aj]))
            || (ctx.md.shake_solute && ai+1 <= ctx.n_atoms_solute) 
            || (ctx.md.shake_solvent && ai+1 > ctx.n_atoms_solute) ) {
            ctx.shake_bonds[shake].ai = ai+1;
            ctx.shake_bonds[shake].aj = aj+1;
            ctx.shake_bonds[shake].dist2 = pow(ctx.cbonds[ctx.bonds[bi].code-1].b0, 2);
            shake++;
        }
    }

    // Get total number of shake constraints in solute (used for separate scaling of temperatures)
    for (int i = 0; i < ctx.n_molecules - ctx.n_waters; i++) {
        n_solute_shake_constraints += ctx.mol_n_shakes[i];
    }

    ctx.Ndegf = 3 * ctx.n_atoms - ctx.n_shake_constraints;
    ctx.Ndegfree = ctx.Ndegf - 3 * ctx.n_excluded + excl_shake;

    ctx.Ndegf_solvent = ctx.Ndegf - 3 * ctx.n_atoms_solute  + n_solute_shake_constraints;
    ctx.Ndegf_solute = ctx.Ndegf - ctx.Ndegf_solvent;

    ctx.Ndegfree_solvent = ctx.Ndegfree - (ctx.n_shake_constraints - n_solute_shake_constraints);
    ctx.Ndegfree_solute = ctx.Ndegfree - ctx.Ndegfree_solvent;

    printf("n_shake_constrains = %d, n_solute_shake_constraints = %d, excl_shake = %f\n", ctx.n_shake_constraints, n_solute_shake_constraints, excl_shake);

    if (ctx.Ndegfree_solvent * ctx.Ndegfree_solute == 0) {
        ctx.separate_scaling = false;
    }
    else {
        ctx.separate_scaling = true;
    }
}

/* =============================================
 * == CALCUTED IN THE INTEGRATION
 * =============================================
 */

void init_patoms() {
    auto &ctx = Context::instance();
    ctx.n_patoms = ctx.n_atoms_solute - ctx.n_qatoms;

    ctx.p_atoms.resize(ctx.n_patoms);

    // Loop through all solutes, adding a p atom to the list every time a non-q atom is encountered
    int pi = 0;
    int qi = 0;
    for (int i = 0; i < ctx.n_atoms_solute; i++) {
        if (ctx.n_qatoms > 0 && qi < ctx.n_qatoms && i == ctx.q_atoms[qi].a-1) {
            qi++;
        }
        else {
            ctx.p_atoms[pi].a = i+1;
            pi++;
        }
    }
}

/* =============================================
 * == ENERGY & TEMPERATURE
 * =============================================
 */


void calc_temperature() {
    auto &ctx = Context::instance();
    printf("Ndegf = %f, Ndegfree = %f, n_excluded = %d, Ndegfree_solvent = %f, Ndegfree_solute = %f\n",
        ctx.Ndegf, ctx.Ndegfree, ctx.n_excluded, ctx.Ndegfree_solvent, ctx.Ndegfree_solute);
    ctx.Temp = 0;
    ctx.Tfree = 0;
    double Temp_solute = 0, Tfree_solute = 0, Texcl_solute = 0;
    double Tfree_solvent = 0, Temp_solvent = 0, Texcl_solvent = 0;
    double Ekinmax = 1000.0 * ctx.Ndegf * Boltz * ctx.md.temperature / 2.0 / ctx.n_atoms;
    double ener;
    double mass_i;

    ctx.Temp = 0;
    for (int i = 0; i < ctx.n_atoms_solute; i++) {
        mass_i = ctx.catypes[ctx.atypes[i].code - 1].m;
        ener = .5 * mass_i * (pow(ctx.velocities[i].x, 2) + pow(ctx.velocities[i].y, 2) + pow(ctx.velocities[i].z, 2));
        Temp_solute += ener;
        if (!ctx.excluded[i]) {
            Tfree_solute += ener;
        }
        else {
            Texcl_solute += ener;
        }
        if (ener > Ekinmax) {
            printf(">>> WARNING: hot atom %d: %f\n", i, ener/Boltz/3);
        }
    }

    for (int i = ctx.n_atoms_solute; i < ctx.n_atoms; i++) {
        mass_i = ctx.catypes[ctx.atypes[i].code - 1].m;
        ener = .5 * mass_i * (pow(ctx.velocities[i].x, 2) + pow(ctx.velocities[i].y, 2) + pow(ctx.velocities[i].z, 2));
        Temp_solvent += ener;
        if (!ctx.excluded[i]) {
            Tfree_solvent += ener;
        }
        else {
            Texcl_solvent += ener;
        }
        if (ener > Ekinmax) {
            printf(">>> WARNING: hot atom %d: %f\n", i, ener/Boltz/3);
        }
    }

    ctx.Tfree = Tfree_solute + Tfree_solvent;
    ctx.Temp = Temp_solute + Temp_solvent;

    ctx.E_total.Ukin = ctx.Temp;

    ctx.Temp = 2.0 * ctx.Temp / Boltz / ctx.Ndegf;
    ctx.Tfree = 2.0 * ctx.Tfree / Boltz / ctx.Ndegfree;

    if (ctx.separate_scaling) {
        if (Tfree_solvent != 0) ctx.Tscale_solvent = sqrt(1 + (ctx.dt / ctx.tau_T) * (ctx.md.temperature / Tfree_solvent - 1.0));
        if (Tfree_solute != 0) ctx.Tscale_solute = sqrt(1 + (ctx.dt / ctx.tau_T) * (ctx.md.temperature / Tfree_solute - 1.0));
    }
    else {
        if (ctx.Tfree != 0) ctx.Tscale_solvent = sqrt(1 + (ctx.dt / ctx.tau_T) * (ctx.md.temperature / ctx.Tfree - 1.0));
        ctx.Tscale_solute = ctx.Tscale_solvent;
    }
    printf("Tscale = %f, tau_T = %f, Temp = %f, Tfree = %f\n", ctx.Tscale_solvent, ctx.tau_T, ctx.Temp, ctx.Tfree);
}

/* =============================================
 * == INTEGRATION METHODS
 * =============================================
 */

void calc_leapfrog() {
    auto &ctx = Context::instance();
    double mass_i, winv_i;
    for (int i = 0; i < ctx.n_atoms_solute; i++) {
        mass_i = ctx.catypes[ctx.atypes[i].code - 1].m;

        winv_i = 1/mass_i;
        ctx.velocities[i].x = (ctx.velocities[i].x - ctx.dvelocities[i].x * ctx.dt * winv_i) * ctx.Tscale_solute;
        ctx.velocities[i].y = (ctx.velocities[i].y - ctx.dvelocities[i].y * ctx.dt * winv_i) * ctx.Tscale_solute;
        ctx.velocities[i].z = (ctx.velocities[i].z - ctx.dvelocities[i].z * ctx.dt * winv_i) * ctx.Tscale_solute;

        // Prepare copy for shake
        ctx.xcoords[i].x = ctx.coords[i].x;
        ctx.xcoords[i].y = ctx.coords[i].y;
        ctx.xcoords[i].z = ctx.coords[i].z;

        ctx.coords[i].x += ctx.velocities[i].x * ctx.dt;
        ctx.coords[i].y += ctx.velocities[i].y * ctx.dt;
        ctx.coords[i].z += ctx.velocities[i].z * ctx.dt;

    }

    for (int i = ctx.n_atoms_solute; i < ctx.n_atoms; i++) {
        mass_i = ctx.catypes[ctx.atypes[i].code - 1].m;

        winv_i = 1/mass_i;
        ctx.velocities[i].x = (ctx.velocities[i].x - ctx.dvelocities[i].x * ctx.dt * winv_i) * ctx.Tscale_solvent;
        ctx.velocities[i].y = (ctx.velocities[i].y - ctx.dvelocities[i].y * ctx.dt * winv_i) * ctx.Tscale_solvent;
        ctx.velocities[i].z = (ctx.velocities[i].z - ctx.dvelocities[i].z * ctx.dt * winv_i) * ctx.Tscale_solvent;

        // Prepare copy for shake
        ctx.xcoords[i].x = ctx.coords[i].x;
        ctx.xcoords[i].y = ctx.coords[i].y;
        ctx.xcoords[i].z = ctx.coords[i].z;

        ctx.coords[i].x += ctx.velocities[i].x * ctx.dt;
        ctx.coords[i].y += ctx.velocities[i].y * ctx.dt;
        ctx.coords[i].z += ctx.velocities[i].z * ctx.dt;

    }


    // Shake if necessary
    if (ctx.n_shake_constraints > 0) {
        calc_shake_constraints(ctx.coords.data(), ctx.xcoords.data());
        for (int i = 0; i < ctx.n_atoms; i++) {
            ctx.velocities[i].x = (ctx.coords[i].x - ctx.xcoords[i].x) / ctx.dt;
            ctx.velocities[i].y = (ctx.coords[i].y - ctx.xcoords[i].y) / ctx.dt;
            ctx.velocities[i].z = (ctx.coords[i].z - ctx.xcoords[i].z) / ctx.dt;
        }
    }
}

// Write header (number of atoms) to output file
void write_header(const char *filename) {
    auto &ctx = Context::instance();
    FILE * fp;

    char path[1024];
    sprintf(path, "%s/output/%s", ctx.base_folder.c_str(), filename);

    fp = fopen(path, "w");
  
    fprintf(fp, "%d\n", ctx.n_atoms);
  
    fclose (fp);
}

// Write header (number of atoms & lambdas) to output file
void write_energy_header() {
    auto &ctx = Context::instance();
    FILE * fp;

    char path[1024];
    sprintf(path, "%s/output/%s", ctx.base_folder.c_str(), "energies.csv");

    fp = fopen(path, "w");
  
    fprintf(fp, "%d\n", ctx.n_atoms);

    fprintf(fp, "lambdas\n");
    fprintf(fp, "%d\n", ctx.n_lambdas);

    for (int state = 0; state < ctx.n_lambdas; state++) {
        fprintf(fp, "%f\n", ctx.lambdas[state]);
    }
  
    fclose (fp);
}

// Write step number, coordinates of atoms to coordinate output file
void write_coords(int iteration) {
    auto &ctx = Context::instance();
    if (iteration % ctx.md.trajectory != 0) return;
    FILE * fp;
    int i;

    char path[1024];
    sprintf(path, "%s/output/%s", ctx.base_folder.c_str(), "coords.csv");

    fp = fopen(path, "a");

    fprintf(fp, "%d\n", iteration / ctx.md.trajectory);
    for(i = 0; i < ctx.n_atoms; i++) {
        fprintf(fp, "%f;%f;%f\n", ctx.coords[i].x, ctx.coords[i].y, ctx.coords[i].z);
    }
  
    fclose (fp);
}

// Write step number, velocities of atoms to coordinate output file
void write_velocities(int iteration) {
    auto &ctx = Context::instance();
    if (iteration % ctx.md.trajectory != 0) return;
    FILE * fp;
    int i;

    char path[1024];
    sprintf(path, "%s/output/%s", ctx.base_folder.c_str(), "velocities.csv");

    fp = fopen(path, "a");

    fprintf(fp, "%d\n", iteration / ctx.md.trajectory);
    for(i = 0; i < ctx.n_atoms; i++) {
        fprintf(fp, "%f;%f;%f\n", ctx.velocities[i].x, ctx.velocities[i].y, ctx.velocities[i].z);
    }
  
    fclose (fp);
}

// Write step number, energies of atoms to coordinate output file
void write_energies(int iteration) {
    auto &ctx = Context::instance();
    if (iteration % ctx.md.energy != 0) return;
    FILE * fp;

    char path[1024];
    sprintf(path, "%s/output/%s", ctx.base_folder.c_str(), "energies.csv");

    fp = fopen(path, "a");

    fprintf(fp, "interval %d\n", iteration / ctx.md.energy);

    fprintf(fp, "[temperature]\n");
    fprintf(fp, "Temp\t%f\n", ctx.Temp);
    fprintf(fp, "\n");

    fprintf(fp, "[bonded]\n");
    fprintf(fp, "type\tUbond\tUangle\tUtor\tUimp\n");
    fprintf(fp, "p\t%f\t%f\t%f\t%f\n", ctx.E_bond_p.Ubond, ctx.E_bond_p.Uangle, ctx.E_bond_p.Utor, ctx.E_bond_p.Uimp);
    fprintf(fp, "w\t%f\t%f\t%f\t%f\n", ctx.E_bond_w.Ubond, ctx.E_bond_w.Uangle, ctx.E_bond_w.Utor, ctx.E_bond_w.Uimp);
    fprintf(fp, "qp\t%f\t%f\t%f\t%f\n", ctx.E_bond_q.Ubond, ctx.E_bond_q.Uangle, ctx.E_bond_q.Utor, ctx.E_bond_q.Uimp);
    fprintf(fp, "\n");

    fprintf(fp, "[nonbonded]\n");
    fprintf(fp, "type\tUcoul\tUvdw\n");
    fprintf(fp, "pp\t%f\t%f\n", ctx.E_nonbond_pp.Ucoul, ctx.E_nonbond_pp.Uvdw);
    fprintf(fp, "pw\t%f\t%f\n", ctx.E_nonbond_pw.Ucoul, ctx.E_nonbond_pw.Uvdw);
    fprintf(fp, "ww\t%f\t%f\n", ctx.E_nonbond_ww.Ucoul, ctx.E_nonbond_ww.Uvdw);
    fprintf(fp, "qx\t%f\t%f\n", ctx.E_nonbond_qx.Ucoul, ctx.E_nonbond_qx.Uvdw);
    fprintf(fp, "\n");

    fprintf(fp, "[restraint]\n");
    fprintf(fp, "Uradx\t%f\n", ctx.E_restraint.Uradx);
    fprintf(fp, "Upolx\t%f\n", ctx.E_restraint.Upolx);
    fprintf(fp, "Ushell\t%f\n", ctx.E_restraint.Ushell);
    fprintf(fp, "Ufix\t%f\n", ctx.E_restraint.Ufix);
    fprintf(fp, "Upres\t%f\n", ctx.E_restraint.Upres);
    fprintf(fp, "Total\t%f\n", ctx.E_restraint.Urestr);
    fprintf(fp, "\n");

    fprintf(fp, "[q-energies]\n");
    fprintf(fp, "lambda\tSUM\tUbond\tUangle\tUtor\tUimp\tUcoul\tUvdw\tUrestr\n");
    for (int state = 0; state < ctx.n_lambdas; state++) {
        fprintf(fp, "%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\n", ctx.lambdas[state], ctx.EQ_total[state].Utot, ctx.EQ_bond[state].Ubond,
            ctx.EQ_bond[state].Uangle, ctx.EQ_bond[state].Utor, ctx.EQ_bond[state].Uimp, ctx.EQ_nonbond_qx[state].Ucoul, ctx.EQ_nonbond_qx[state].Uvdw, ctx.EQ_restraint[state].Urestr);
    }
    fprintf(fp, "\n");
    
    fprintf(fp, "[total]\n");
    fprintf(fp, "Ukin\t%f\n", ctx.E_total.Ukin);
    fprintf(fp, "Upot\t%f\n", ctx.E_total.Upot);
    fprintf(fp, "Utot\t%f\n", ctx.E_total.Utot);
    fprintf(fp, "\n");
  
    fclose (fp);
}

void cpu_run(int n_steps) {
    for (int i = 0; i <= n_steps; i++) {
        calc_integration_step(i);
    }
}

void calc_integration() {
    auto &ctx = Context::instance();
    init_variables();
    if (ctx.run_gpu) {
        auto &handler = CudaHandler::instance();
        handler.initialize();
        handler.run(ctx.md.steps + 1);
        handler.shutdown();
    } else {
        cpu_run(ctx.md.steps);
    }
    clean_variables();
}

void reset_energies() {
    auto &ctx = Context::instance();
    for (int i = 0; i < ctx.n_atoms; i++) {
        ctx.dvelocities[i].x = 0;
        ctx.dvelocities[i].y = 0;
        ctx.dvelocities[i].z = 0;
    }
    ctx.E_total.Upot = 0;
    ctx.E_bond_p.Uangle = 0;
    ctx.E_bond_p.Ubond = 0;
    ctx.E_bond_p.Utor = 0;
    ctx.E_bond_p.Uimp = 0;
    ctx.E_bond_w.Uangle = 0;
    ctx.E_bond_w.Ubond = 0;
    ctx.E_bond_w.Utor = 0;
    ctx.E_bond_w.Uimp = 0;
    ctx.E_bond_q.Uangle = 0;
    ctx.E_bond_q.Ubond = 0;
    ctx.E_bond_q.Utor = 0;
    ctx.E_bond_q.Uimp = 0;
    ctx.E_nonbond_pp.Ucoul = 0;
    ctx.E_nonbond_pp.Uvdw = 0;
    ctx.E_nonbond_pw.Ucoul = 0;
    ctx.E_nonbond_pw.Uvdw = 0;
    ctx.E_nonbond_ww.Ucoul = 0;
    ctx.E_nonbond_ww.Uvdw = 0;
    ctx.E_nonbond_qx.Ucoul = 0;
    ctx.E_nonbond_qx.Uvdw = 0;
    ctx.E_restraint.Uradx = 0;
    ctx.E_restraint.Upolx = 0;
    ctx.E_restraint.Ufix = 0;
    ctx.E_restraint.Ushell = 0;
    ctx.E_restraint.Upres = 0;
    ctx.E_restraint.Urestr = 0;
    for (int state = 0; state < ctx.n_lambdas; state++) {
        ctx.EQ_bond[state].Uangle = 0;
        ctx.EQ_bond[state].Ubond = 0;
        ctx.EQ_bond[state].Utor = 0;
        ctx.EQ_bond[state].Uimp = 0;
        ctx.EQ_nonbond_qq[state].Ucoul = 0;
        ctx.EQ_nonbond_qq[state].Uvdw = 0;
        ctx.EQ_nonbond_qp[state].Ucoul = 0;
        ctx.EQ_nonbond_qp[state].Uvdw = 0;
        ctx.EQ_nonbond_qw[state].Ucoul = 0;
        ctx.EQ_nonbond_qw[state].Uvdw = 0;
        ctx.EQ_restraint[state].Urestr = 0;
    }
}

void calc_bonded_forces() {
    auto &ctx = Context::instance();
    ctx.E_bond_p.Uangle = calc_angle_forces(0, ctx.n_angles_solute);
    ctx.E_bond_w.Uangle = calc_angle_forces(ctx.n_angles_solute, ctx.n_angles);

    ctx.E_bond_p.Ubond = calc_bond_forces(0, ctx.n_bonds_solute);
    ctx.E_bond_w.Ubond = calc_bond_forces(ctx.n_bonds_solute, ctx.n_bonds);

    ctx.E_bond_p.Utor = calc_torsion_forces(0, ctx.n_torsions_solute);
    ctx.E_bond_w.Utor = calc_torsion_forces(ctx.n_torsions_solute, ctx.n_torsions);

    ctx.E_bond_p.Uimp = calc_improper2_forces(0, ctx.n_impropers_solute);
    ctx.E_bond_w.Uimp = calc_improper2_forces(ctx.n_impropers_solute, ctx.n_impropers);
}

void calc_integration_step(int iteration) {
    auto &ctx = Context::instance();
    printf("================================================\n");
    if (iteration > 0) {
        printf("== STEP %d\n", iteration);
    }
    else {
        printf("== INITIAL ENERGIES\n");
    }
    printf("================================================\n");

    // Reset derivatives & energies
    reset_energies();

    // Determine temperature and kinetic energy
    calc_temperature();

    // Determine acceleration
    clock_t start = clock();

    // First solute interactions
    calc_bonded_forces();

    clock_t end_bonded = clock();

    clock_t start_pp, end_pp, start_qp, end_qp;
    start_qp = clock();
    calc_nonbonded_qp_forces();
    end_qp = clock();
    start_pp = clock();
    calc_nonbonded_pp_forces();
    end_pp = clock();

    clock_t start_ww, end_ww, start_pw, end_pw;
    // Now solvent interactions
    if (ctx.n_waters > 0) {
        start_ww = clock();
        calc_nonbonded_ww_forces();
        end_ww = clock();
        start_pw = clock();
        calc_nonbonded_pw_forces();
        end_pw = clock();
        calc_nonbonded_qw_forces();
    }

    clock_t end_nonbonded = clock();

    // Calculate restraints
    if (ctx.n_waters > 0) {
        calc_radix_w_forces();
        if (ctx.md.polarisation) {
            calc_polx_w_forces(iteration);
        }
    }

    calc_pshell_forces();
    calc_restrseq_forces();
    calc_restrdis_forces();
    calc_restrpos_forces();
    calc_restrang_forces();
    calc_restrwall_forces();
    
    
    clock_t end_restraints = clock();

    // Q-Q nonbonded interactions
    clock_t start_qq = clock();

    calc_nonbonded_qq_forces();
    clock_t end_qq = clock();

    // Q-atom bonded interactions: loop over Q-atom states
    for (int state = 0; state < ctx.n_lambdas; state++) {
        calc_qangle_forces(state);
        calc_qbond_forces(state);
        calc_qtorsion_forces(state);
    }

    clock_t end_qatoms = clock();

    // Now apply leapfrog integration
    calc_leapfrog();

    // Recalculate temperature and kinetic energy for output
    calc_temperature();
    // Update total potential energies with an average of all states
    for (int state = 0; state < ctx.n_lambdas; state++) {
        if (ctx.lambdas[state] == 0) {
            ctx.EQ_bond[state].Uangle = 0;
            ctx.EQ_bond[state].Ubond = 0;
            ctx.EQ_bond[state].Utor = 0;
            ctx.EQ_bond[state].Uimp = 0;
            ctx.EQ_nonbond_qq[state].Ucoul = 0;
            ctx.EQ_nonbond_qq[state].Uvdw = 0;
            ctx.EQ_nonbond_qp[state].Ucoul = 0;
            ctx.EQ_nonbond_qp[state].Uvdw = 0;
            ctx.EQ_nonbond_qw[state].Ucoul = 0;
            ctx.EQ_nonbond_qw[state].Uvdw = 0;
            ctx.EQ_restraint[state].Urestr = 0;
        }

        ctx.EQ_nonbond_qx[state].Ucoul = ctx.EQ_nonbond_qq[state].Ucoul + ctx.EQ_nonbond_qp[state].Ucoul + ctx.EQ_nonbond_qw[state].Ucoul;
        ctx.EQ_nonbond_qx[state].Uvdw = ctx.EQ_nonbond_qq[state].Uvdw + ctx.EQ_nonbond_qp[state].Uvdw + ctx.EQ_nonbond_qw[state].Uvdw;

        ctx.EQ_total[state].Utot = ctx.EQ_bond[state].Ubond + ctx.EQ_bond[state].Uangle + ctx.EQ_bond[state].Utor + ctx.EQ_bond[state].Uimp
            + ctx.EQ_nonbond_qx[state].Ucoul + ctx.EQ_nonbond_qx[state].Uvdw + ctx.EQ_restraint[state].Urestr;

        ctx.E_bond_q.Ubond += ctx.EQ_bond[state].Ubond * ctx.lambdas[state];
        ctx.E_bond_q.Uangle += ctx.EQ_bond[state].Uangle * ctx.lambdas[state];
        ctx.E_bond_q.Utor += ctx.EQ_bond[state].Utor * ctx.lambdas[state];
        ctx.E_bond_q.Uimp += ctx.EQ_bond[state].Uimp * ctx.lambdas[state];
        ctx.E_nonbond_qx.Ucoul += ctx.EQ_nonbond_qx[state].Ucoul * ctx.lambdas[state];
        ctx.E_nonbond_qx.Uvdw += ctx.EQ_nonbond_qx[state].Uvdw * ctx.lambdas[state];

        // Update protein restraint energies with an average of all states
        ctx.E_restraint.Upres += ctx.EQ_restraint[state].Urestr * ctx.lambdas[state];
    }

    // Update totals
    ctx.E_restraint.Urestr = ctx.E_restraint.Uradx + ctx.E_restraint.Upolx + ctx.E_restraint.Ushell + ctx.E_restraint.Ufix + ctx.E_restraint.Upres;
    ctx.E_total.Upot = ctx.E_bond_p.Ubond + ctx.E_bond_w.Ubond + ctx.E_bond_p.Uangle + ctx.E_bond_w.Uangle + ctx.E_bond_p.Utor + ctx.E_bond_p.Uimp
        + ctx.E_nonbond_pp.Ucoul + ctx.E_nonbond_pp.Uvdw + ctx.E_nonbond_pw.Ucoul + ctx.E_nonbond_pw.Uvdw + ctx.E_nonbond_ww.Ucoul
        + ctx.E_nonbond_ww.Uvdw + ctx.E_bond_q.Ubond + ctx.E_bond_q.Uangle + ctx.E_bond_q.Utor + ctx.E_bond_q.Uimp
        + ctx.E_nonbond_qx.Ucoul + ctx.E_nonbond_qx.Uvdw + ctx.E_restraint.Urestr;
    ctx.E_total.Utot = ctx.E_total.Upot + ctx.E_total.Ukin;

    printf("[temperature]\n");
    printf("Temp\t%f\n", ctx.Temp);
    printf("\n");

    printf("[bonded]\n");
    printf("type\tUbond\tUangle\tUtor\tUimp\n");
    printf("p\t%f\t%f\t%f\t%f\n", ctx.E_bond_p.Ubond, ctx.E_bond_p.Uangle, ctx.E_bond_p.Utor, ctx.E_bond_p.Uimp);
    printf("w\t%f\t%f\t%f\t%f\n", ctx.E_bond_w.Ubond, ctx.E_bond_w.Uangle, ctx.E_bond_w.Utor, ctx.E_bond_w.Uimp);
    printf("qp\t%f\t%f\t%f\t%f\n", ctx.E_bond_q.Ubond, ctx.E_bond_q.Uangle, ctx.E_bond_q.Utor, ctx.E_bond_q.Uimp);
    printf("\n");

    printf("[nonbonded]\n");
    printf("type\tUcoul\tUvdw\n");
    printf("pp\t%f\t%f\n", ctx.E_nonbond_pp.Ucoul, ctx.E_nonbond_pp.Uvdw);
    printf("pw\t%f\t%f\n", ctx.E_nonbond_pw.Ucoul, ctx.E_nonbond_pw.Uvdw);
    printf("ww\t%f\t%f\n", ctx.E_nonbond_ww.Ucoul, ctx.E_nonbond_ww.Uvdw);
    printf("qx\t%f\t%f\n", ctx.E_nonbond_qx.Ucoul, ctx.E_nonbond_qx.Uvdw);
    printf("\n");

    printf("[restraint]\n");
    printf("Uradx\t%f\n", ctx.E_restraint.Uradx);
    printf("Upolx\t%f\n", ctx.E_restraint.Upolx);
    printf("Ushell\t%f\n", ctx.E_restraint.Ushell);
    printf("Ufix\t%f\n", ctx.E_restraint.Ufix);
    printf("Upres\t%f\n", ctx.E_restraint.Upres);
    printf("Total\t%f\n", ctx.E_restraint.Urestr);
    printf("\n");

    printf("[q-energies]\n");
    printf("lambda\tSUM\tUbond\tUangle\tUtor\tUimp\tUcoul\tUvdw\tUrestr\n");
    for (int state = 0; state < ctx.n_lambdas; state++) {
        printf("%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\n", ctx.lambdas[state], ctx.EQ_total[state].Utot, ctx.EQ_bond[state].Ubond,
            ctx.EQ_bond[state].Uangle, ctx.EQ_bond[state].Utor, ctx.EQ_bond[state].Uimp, ctx.EQ_nonbond_qx[state].Ucoul, ctx.EQ_nonbond_qx[state].Uvdw, ctx.EQ_restraint[state].Urestr);
    }
    printf("\n");
    
    printf("[total]\n");
    printf("Ukin\t%f\n", ctx.E_total.Ukin);
    printf("Upot\t%f\n", ctx.E_total.Upot);
    printf("Utot\t%f\n", ctx.E_total.Utot);
    printf("\n");

    // Append output files
    write_coords(iteration);
    write_velocities(iteration);
    write_energies(iteration);    

    clock_t end_calculation = clock();

    // Profiler info
#ifdef __PROFILING__
    printf("Elapsed time for bonded forces: %f\n", (end_bonded-start) / (double)CLOCKS_PER_SEC );
    printf("Elapsed time for non-bonded forces: %f\n", (end_nonbonded-end_bonded) / (double)CLOCKS_PER_SEC);
    printf("Elapsed time for pp interactions: %f\n", (end_pp-start_pp) / (double)CLOCKS_PER_SEC );
    printf("Elapsed time for qq interaction: %f\n",  (end_qq-start_qq) / (double)CLOCKS_PER_SEC );
    printf("Elapsed time for qp interaction: %f\n",  (end_qp-start_qp) / (double)CLOCKS_PER_SEC );
    if (ctx.n_waters > 0) {
        printf("Elapsed time for ww interactions: %f\n", (end_ww-start_ww) / (double)CLOCKS_PER_SEC );
        printf("Elapsed time for pw interactions: %f\n", (end_pw-start_pw) / (double)CLOCKS_PER_SEC );
    }
    printf("---\n");
    printf("Elapsed time for entire time-step: %f\n", (end_calculation-start) / (double)CLOCKS_PER_SEC);
#endif /* __PROFILING__ */

}


void init_variables() {
    auto &ctx = Context::instance();
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
    //exclude_all_atoms_excluded_definitions();
    
    // Shake constraints, need to be initialized before last part of shrink_topology
    if (ctx.md.charge_groups) {
        if (ctx.iuse_switch_atom == 1) {
            init_pshells_with_switch_atoms();
        } else {
            init_pshells_with_centroids();
        }
    }
    else {
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
    auto &ctx = Context::instance();

    for (auto &charge_group : ctx.charge_groups) {
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
