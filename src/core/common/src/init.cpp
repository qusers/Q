#include <algorithm>
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

static void finalize_ngbrs14_host() {
    auto& ctx = Context::instance();
    auto &LJ_matrix = ctx.LJ_matrix->cpu_data_p;
    std::vector<int3> ngbrs_14;
    ngbrs_14.reserve(ctx.n_atoms_solute);

    for (int i = 0; i < ctx.n_atoms_solute; i++) {
        for (int j = i + 1; j < ctx.n_atoms_solute; j++) {
            if (LJ_matrix[i * ctx.n_atoms_solute + j] == 1) {
                int pair_type = NONBONDED_14_PP;
                const bool ai_is_q = ctx.atom_to_qi[i] != -1;
                const bool aj_is_q = ctx.atom_to_qi[j] != -1;
                if (ai_is_q && aj_is_q) {
                    pair_type = NONBONDED_14_QQ;
                } else if (ai_is_q || aj_is_q) {
                    pair_type = NONBONDED_14_QP;
                }
                ngbrs_14.push_back({i, j, pair_type});
            }
        }
    }
    ctx.ngbrs_14 = std::make_unique<HostDeviceBuffer<int3>>(ngbrs_14.size(), true, ctx.run_gpu);
    if (!ngbrs_14.empty()) {
        std::copy(ngbrs_14.begin(), ngbrs_14.end(), ctx.ngbrs_14->cpu_data_p);
    }
    if (ctx.run_gpu) {
        ctx.LJ_matrix->upload();
        ctx.ngbrs_14->upload();
    }
    ctx.n_ngbrs14 = static_cast<int>(ctx.ngbrs_14->length);
}

// Remove bonds, angles, torsions and impropers which are excluded or changed in the FEP file
void exclude_qatom_definitions() {
    auto& ctx = Context::instance();
    int excluded;
    int ai = 0, bi = 0, ii = 0, ti = 0;
    int qai = 0, qbi = 0, qii = 0, qti = 0;

    excluded = 0;
    int solute_excluded = 0;

    auto &bonds = ctx.bonds->cpu_data_p;
    auto &angles = ctx.angles->cpu_data_p;
    auto &torsions = ctx.torsions->cpu_data_p;

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
    auto *excluded = ctx.excluded->cpu_data_p;
    int n_excl;
    int ai = 0, bi = 0, ii = 0, ti = 0;
    auto &impropers = ctx.impropers->cpu_data_p;
    auto &torsions = ctx.torsions->cpu_data_p;

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

void exclude_shaken_definitions() {
    auto& ctx = Context::instance();
    int excluded;
    int solute_excluded;
    int bi = 0;
    int si = 0;
    int ang_i = 0;
    int ai, aj;
    auto &shake_bonds = ctx.shake_bonds->cpu_data_p;

    excluded = 0;
    solute_excluded = 0;
    auto &bonds = ctx.bonds->cpu_data_p;
    if (ctx.n_shake_constraints > 0) {
        for (int i = 0; i < ctx.n_bonds; i++) {
            if (bonds[i].ai == shake_bonds[si].ai && bonds[i].aj == shake_bonds[si].aj) {
                si++;
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

    excluded = 0;
    solute_excluded = 0;
    auto &angles = ctx.angles->cpu_data_p;
    if (ctx.n_shake_constraints > 0) {
        // Mark angles whose terminal atoms (i and k) match a shaken bond
        for (int i = 0; i < ctx.n_shake_constraints; i++) {
            ai = shake_bonds[i].ai;
            aj = shake_bonds[i].aj;
            for (int j = 0; j < ctx.n_angles; j++) {
                if ((angles[j].ai == ai && angles[j].ak == aj) || (angles[j].ai == aj && angles[j].ak == ai)) {
                    angles[j].code = 0;
                    break;
                }
            }
        }

        for (int i = 0; i < ctx.n_angles; i++) {
            if (angles[i].code == 0) {
                excluded++;
                if (i < ctx.n_angles_solute) solute_excluded++;
            } else {
                angles[ang_i] = angles[i];
                ang_i++;
            }
        }
        ctx.n_angles -= excluded;
        ctx.n_angles_solute -= solute_excluded;
    }
}

void init_velocities() {
    auto& ctx = Context::instance();
    auto &atypes = ctx.atypes->cpu_data_p;
    auto &catypes = ctx.catypes->cpu_data_p;
    ctx.velocities = std::make_unique<HostDeviceBuffer<vel_t>>(ctx.n_atoms, true, ctx.run_gpu);
    auto &velocities = ctx.velocities->cpu_data_p;

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

    if (ctx.run_gpu) {
        ctx.velocities->upload();
    }
}

void init_dvelocities() {
    auto& ctx = Context::instance();
    ctx.dvelocities = std::make_unique<HostDeviceBuffer<dvel_t>>(ctx.n_atoms, true, ctx.run_gpu);
}

void init_xcoords() {
    auto& ctx = Context::instance();
    ctx.xcoords = std::make_unique<HostDeviceBuffer<coord_t>>(ctx.n_atoms, true, ctx.run_gpu);
}

void init_inv_mass() {
    auto& ctx = Context::instance();
    auto &atypes = ctx.atypes->cpu_data_p;
    auto &catypes = ctx.catypes->cpu_data_p;
    ctx.winv = std::make_unique<HostDeviceBuffer<double>>(ctx.n_atoms, true, ctx.run_gpu);
    auto *winv = ctx.winv->cpu_data_p;
    for (int ai = 0; ai < ctx.n_atoms; ai++) {
        winv[ai] = 1 / catypes[atypes[ai].code - 1].m;
    }

    if (ctx.run_gpu) {
        ctx.winv->upload();
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
    auto &bonds = ctx.bonds->cpu_data_p;
    auto &cbonds = ctx.cbonds->cpu_data_p;
    auto &angles = ctx.angles->cpu_data_p;
    auto &cangles = ctx.cangles->cpu_data_p;
    // Get water properties from the first water molecule.
    cbond_t cbondw = cbonds[bonds[ctx.n_atoms_solute].code - 1];
    cangle_t canglew = cangles[angles[ctx.n_atoms_solute].code - 1];
    const double crg_ow = ctx.unified_ccharge(ctx.n_atoms_solute, 0).charge;
    const double mu_w = -crg_ow * cbondw.b0 * cos(canglew.th0 / 2);

    drs = wpolr_layer / drouter;

    ctx.n_shells = (int)floor(-0.5 + sqrt(2 * drs + 0.25));
    ctx.wshells = std::make_unique<HostDeviceBuffer<shell_t>>(ctx.n_shells, true, ctx.run_gpu);
    auto *wshells = ctx.wshells->cpu_data_p;

    printf("n_shells = %d\n", ctx.n_shells);

    router = ctx.topo.solvent_radius;
    ctx.n_max_inshell = 0;

    for (int i = 0; i < ctx.n_shells; i++) {
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

    // rc > wshells[n_shells-1].router - wshells[n_shells-1].dr
    printf("shell 0: (%f, %f). shell 1: (%f, %f). shell 2: (%f, %f).\n", wshells[0].router, wshells[0].router - wshells[0].dr, wshells[1].router, wshells[1].router - wshells[1].dr, wshells[2].router, wshells[2].router - wshells[2].dr);

    ctx.n_max_inshell = ctx.n_waters;  // Make largest a little bigger just in case

    // Initialize arrays needed for bookkeeping
    ctx.theta.assign(ctx.n_waters, 0.0);
    ctx.theta0.assign(ctx.n_waters, 0.0);
    ctx.tdum.assign(ctx.n_waters, 0.0);

    ctx.list_sh.assign(ctx.n_max_inshell, std::vector<int>(ctx.n_shells, 0));
    ctx.nsort.assign(ctx.n_max_inshell, std::vector<int>(ctx.n_shells, 0));

    if (ctx.run_gpu) {
        ctx.wshells->upload();
    }
}

void init_pshells() {
    auto& ctx = Context::instance();
    auto &atypes = ctx.atypes->cpu_data_p;
    auto &catypes = ctx.catypes->cpu_data_p;
    auto &coords_init = ctx.coords_init->cpu_data_p;
    auto *excluded = ctx.excluded->cpu_data_p;
    double mass, r2, rin2;

    ctx.heavy = std::make_unique<HostDeviceBuffer<bool>>(ctx.n_atoms, true, ctx.run_gpu);
    auto *heavy = ctx.heavy->cpu_data_p;
    ctx.shell = std::make_unique<HostDeviceBuffer<bool>>(ctx.n_atoms, true, ctx.run_gpu);
    auto *shell = ctx.shell->cpu_data_p;
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

    if (ctx.run_gpu) {
        ctx.heavy->upload();
        ctx.shell->upload();
    }

    printf("n_heavy = %d, n_inshell = %d\n", n_heavy, n_inshell);
}

// Marks heavy atoms for the shell/excluded arrays.
// Shared between switch-atom and centroid shell init paths.
static int mark_heavy_atoms(Context& ctx) {
    auto &atypes = ctx.atypes->cpu_data_p;
    auto &catypes = ctx.catypes->cpu_data_p;
    auto *heavy = ctx.heavy->cpu_data_p;
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
    ctx.heavy = std::make_unique<HostDeviceBuffer<bool>>(ctx.n_atoms, true, ctx.run_gpu);
    auto *heavy = ctx.heavy->cpu_data_p;
    ctx.shell = std::make_unique<HostDeviceBuffer<bool>>(ctx.n_atoms, true, ctx.run_gpu);
    auto *shell = ctx.shell->cpu_data_p;
    for (int i = 0; i < ctx.n_atoms; ++i) {
        heavy[i] = false;
        shell[i] = false;
    }
}

static void init_pshells_from_charge_groups(const charge_group_config_t& charge_groups) {
    auto& ctx = Context::instance();
    auto &coords_init = ctx.coords_init->cpu_data_p;
    auto *excluded = ctx.excluded->cpu_data_p;
    double r2, rin2;
    const bool use_switch_atom = charge_groups.iuse_switch_atom == 1;

    allocate_pshell_buffers(ctx);
    auto *heavy = ctx.heavy->cpu_data_p;
    auto *shell = ctx.shell->cpu_data_p;
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

    if (ctx.run_gpu) {
        ctx.heavy->upload();
        ctx.shell->upload();
    }

    printf("(%s): n_heavy = %d, n_inshell = %d\n", use_switch_atom ? "switch atoms" : "centroids", n_heavy, n_inshell);
}

void init_restrseqs() {
    auto& ctx = Context::instance();
    ctx.n_restrseqs = 1;
    ctx.restrseqs = std::make_unique<HostDeviceBuffer<restrseq_t>>(1, true, ctx.run_gpu);

    restrseq_t seq;
    seq.ai = 1;
    seq.aj = 14;
    seq.k = 1.0;
    seq.ih = 0;
    seq.to_center = 2;

    ctx.restrseqs->cpu_data_p[0] = seq;

    if (ctx.run_gpu) {
        ctx.restrseqs->upload();
    }
}

/* =============================================
 * == SHAKE
 * =============================================
 */

void init_shake() {
    auto& ctx = Context::instance();
    auto *heavy = ctx.heavy->cpu_data_p;
    auto *excluded = ctx.excluded->cpu_data_p;
    int ai, aj;
    int mol = 0;
    int shake;
    int n_solute_shake_constraints = 0;
    double excl_shake = 0;
    auto &bonds = ctx.bonds->cpu_data_p;
    auto &cbonds = ctx.cbonds->cpu_data_p;

    ctx.n_shake_constraints = 0;
    ctx.mol_n_shakes = std::make_unique<HostDeviceBuffer<int>>(ctx.n_molecules, true, ctx.run_gpu);
    auto *mol_n_shakes = ctx.mol_n_shakes->cpu_data_p;

    for (int bi = 0; bi < ctx.n_bonds; bi++) {
        ai = bonds[bi].ai - 1;
        aj = bonds[bi].aj - 1;

        while (mol + 1 < ctx.n_molecules && ai + 1 >= ctx.molecules[mol + 1]) {
            // new molecule
            mol += 1;
        }

        if ((ctx.md.shake_hydrogens && (!heavy[ai] || !heavy[aj])) || (ctx.md.shake_solute && ai + 1 <= ctx.n_atoms_solute) || (ctx.md.shake_solvent && ai + 1 > ctx.n_atoms_solute)) {
            mol_n_shakes[mol]++;
            ctx.n_shake_constraints++;

            if (excluded[ai]) excl_shake += 0.5;
            if (excluded[aj]) excl_shake += 0.5;
        }
    }

    ctx.shake_bonds = std::make_unique<HostDeviceBuffer<shake_bond_t>>(ctx.n_shake_constraints, true, ctx.run_gpu);
    auto *shake_bonds = ctx.shake_bonds->cpu_data_p;
    mol = 0;
    shake = 0;
    for (int bi = 0; bi < ctx.n_bonds; bi++) {
        ai = bonds[bi].ai - 1;
        aj = bonds[bi].aj - 1;

        while (mol + 1 < ctx.n_molecules && ai + 1 >= ctx.molecules[mol + 1]) {
            // new molecule
            mol += 1;
        }

        if ((ctx.md.shake_hydrogens && (!heavy[ai] || !heavy[aj])) || (ctx.md.shake_solute && ai + 1 <= ctx.n_atoms_solute) || (ctx.md.shake_solvent && ai + 1 > ctx.n_atoms_solute)) {
            shake_bonds[shake].ai = ai + 1;
            shake_bonds[shake].aj = aj + 1;
            shake_bonds[shake].dist2 = pow(cbonds[bonds[bi].code - 1].b0, 2);
            shake++;
        }
    }

    // Get total number of shake constraints in solute (used for separate scaling of temperatures)
    for (int i = 0; i < ctx.n_molecules - ctx.n_waters; i++) {
        n_solute_shake_constraints += mol_n_shakes[i];
    }

    if (ctx.run_gpu) {
        ctx.mol_n_shakes->upload();
        ctx.shake_bonds->upload();
    }

    ctx.Ndegf = 3 * ctx.n_atoms - ctx.n_shake_constraints;
    ctx.Ndegfree = ctx.Ndegf - 3 * ctx.n_excluded + excl_shake;

    const double Ndegf_solvent = ctx.Ndegf - 3 * ctx.n_atoms_solute + n_solute_shake_constraints;

    const double Ndegfree_solvent = ctx.Ndegfree - (ctx.n_shake_constraints - n_solute_shake_constraints);
    const double Ndegfree_solute = ctx.Ndegfree - Ndegfree_solvent;

    printf("n_shake_constrains = %d, n_solute_shake_constraints = %d, excl_shake = %f\n", ctx.n_shake_constraints, n_solute_shake_constraints, excl_shake);

    if (Ndegfree_solvent * Ndegfree_solute == 0) {
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
static void init_unified_atom_parameters() {
    auto& ctx = Context::instance();
    auto &charges = ctx.charges->cpu_data_p;
    auto &ccharges = ctx.ccharges->cpu_data_p;
    auto &atypes = ctx.atypes->cpu_data_p;
    auto &catypes = ctx.catypes->cpu_data_p;
    const int n_states = ctx.n_parameter_states();

    ctx.atom_to_qi.assign(ctx.n_atoms, -1);
    for (int qi = 0; qi < ctx.n_qatoms; qi++) {
        ctx.atom_to_qi[ctx.q_atoms[qi]] = qi;
    }

    const int n_extra_q_states = n_states > 0 ? n_states - 1 : 0;
    const int n_codes = ctx.n_atoms + ctx.n_qatoms * n_extra_q_states;
    ctx.unified_ccharges = std::make_unique<HostDeviceBuffer<ccharge_t>>(n_codes, true, ctx.run_gpu);
    ctx.unified_catypes = std::make_unique<HostDeviceBuffer<catype_t>>(n_codes, true, ctx.run_gpu);
    auto *unified_ccharges = ctx.unified_ccharges->cpu_data_p;
    auto *unified_catypes = ctx.unified_catypes->cpu_data_p;

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

    if (ctx.run_gpu) {
        ctx.unified_ccharges->upload();
        ctx.unified_catypes->upload();
    }
}

void init_variables() {
    auto& ctx = Context::instance();
    // From MD file
    parse_md("md.csv");

    ctx.dt = time_unit * ctx.md.stepsize;
    ctx.tau_T = time_unit * ctx.md.bath_coupling;

    if (ctx.run_gpu && ctx.n_lambdas > 2) {
        printf(">>> FATAL: More than 2 states not supported on GPU architecture. Exiting...\n");
        exit(EXIT_FAILURE);
    }

    // From topology file
    parse_topo("topo.csv");

    parse_angles("angles.csv");
    parse_atypes("atypes.csv");
    parse_bonds("bonds.csv");
    parse_cangles("cangles.csv");
    parse_catypes("catypes.csv");
    parse_cbonds("cbonds.csv");
    parse_ccharges("ccharges.csv");
    parse_charges("charges.csv");
    parse_cimpropers("cimpropers.csv");
    parse_coords("coords.csv");
    parse_ctorsions("ctorsions.csv");
    parse_excluded("excluded.csv");
    parse_molecules("molecules.csv");
    parse_impropers("impropers.csv");
    parse_torsions("torsions.csv");
    parse_LJ_matrix();
    parse_ngbrs14("ngbrs14.csv");
    parse_ngbrs23("ngbrs23.csv");
    parse_ngbrs14_long("ngbrs14long.csv");
    parse_ngbrs23_long("ngbrs23long.csv");
    // init_restrseqs();
    init_inv_mass();
    // From FEP file
    parse_qatoms("q_atoms.csv");
    parse_qangcouples("q_angcouples.csv");
    parse_qcangles("q_cangles.csv");
    parse_qcatypes("q_catypes.csv");
    parse_qcbonds("q_cbonds.csv");
    parse_qcimpropers("q_cimpropers.csv");
    parse_qctorsions("q_ctorsions.csv");
    parse_qoffdiags("q_offdiags.csv");
    parse_qimprcouples("q_imprcouples.csv");
    parse_qsoftpairs("q_softpairs.csv");
    parse_qtorcouples("q_torcouples.csv");

    parse_qangles("q_angles.csv");
    parse_qatypes("q_atypes.csv");
    parse_qbonds("q_bonds.csv");
    parse_qcharges("q_charges.csv");
    parse_qelscales("q_elscales.csv");
    parse_qexclpairs("q_exclpairs.csv");
    parse_qimpropers("q_impropers.csv");
    parse_qshakes("q_shakes.csv");
    parse_qsoftcores("q_softcores.csv");
    parse_qtorsions("q_torsions.csv");

    // First part of shrink topology, this needs to be done first as shake constraints are based on bonds
    exclude_qatom_definitions();
    // exclude_all_atoms_excluded_definitions();

    // Shake constraints, need to be initialized before last part of shrink_topology
    if (ctx.md.charge_groups) {
        init_pshells_from_charge_groups(read_charge_groups("charge_groups.csv"));
    } else {
        init_pshells();
    }
    init_shake();

    // Now remove shaken bonds
    exclude_shaken_definitions();

    init_unified_atom_parameters();

    // Init random seed from MD file
    srand(ctx.md.random_seed);

    // From calculation in the integration
    init_patoms();
    init_velocities();
    init_dvelocities();
    init_xcoords();

    // From input file
    parse_icoords("i_coords.csv");
    parse_ivelocities("i_velocities.csv");

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
    ctx.EQ_restraint = std::make_unique<HostDeviceBuffer<E_restraint_t>>(ctx.n_lambdas, true, ctx.run_gpu);

    if (ctx.n_shake_constraints > 0) {
        initial_shaking();
        stop_cm_translation();
    }

    // Write header to file
    write_header("coords.csv");
    write_header("velocities.csv");
    write_energy_header();

    finalize_ngbrs14_host();
}

void clean_variables() {
    auto& ctx = Context::instance();
    ctx.angles.reset();
    ctx.cangles.reset();
    ctx.bonds.reset();
    ctx.cbonds.reset();
    ctx.coords_init.reset();

    ctx.atypes.reset();
    ctx.catypes.reset();
    ctx.ccharges.reset();
    ctx.charges.reset();
    ctx.atom_to_qi.clear();
    ctx.unified_ccharges.reset();
    ctx.unified_catypes.reset();
    ctx.ngbrs_14.reset();
    ctx.p_atoms_list.reset();
    ctx.w_atoms_list.reset();
    ctx.q_atoms_list.reset();
    ctx.charge_pair_products.reset();
    ctx.p_charge_types.reset();
    ctx.w_charge_types.reset();
    ctx.q_charge_types.reset();
    ctx.catype_pair_params.reset();
    ctx.p_catype_types.reset();
    ctx.w_catype_types.reset();
    ctx.q_catype_types.reset();
    ctx.n_charge_types = 0;
    ctx.zero_charge_type = -1;
    ctx.n_catype_types = 0;
    ctx.zero_catype_type = -1;
    ctx.cimpropers.reset();
    ctx.coords.reset();
    ctx.ctorsions.reset();
    ctx.excluded.reset();
    ctx.heavy.reset();
    ctx.impropers.reset();
    ctx.winv.reset();
    ctx.torsions.reset();
    ctx.LJ_matrix.reset();
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
    ctx.q_elscales.reset();
    ctx.q_exclpairs.clear();
    ctx.q_impropers.clear();
    ctx.q_shakes.clear();
    ctx.q_softcores.clear();
    ctx.q_torsions.clear();
    ctx.lambdas.reset();
    ctx.wshells.reset();
    ctx.theta.clear();
    ctx.theta0.clear();
    ctx.tdum.clear();
    ctx.list_sh.clear();
    ctx.nsort.clear();
    ctx.restrseqs.reset();
    ctx.restrangs.reset();
    ctx.restrdists.reset();
    ctx.restrspos.reset();
    ctx.restrwalls.reset();
    ctx.shell.reset();
    ctx.velocities.reset();
    ctx.dvelocities.reset();
    ctx.xcoords.reset();
    ctx.mol_n_shakes.reset();
    ctx.shake_bonds.reset();
    ctx.EQ_total.clear();
    ctx.EQ_bond.clear();
    ctx.EQ_nonbond_qq.clear();
    ctx.EQ_nonbond_qp.clear();
    ctx.EQ_nonbond_qw.clear();
    ctx.EQ_nonbond_qx.clear();
    ctx.EQ_restraint.reset();
}
