#include "cpu_shake.h"

#include <cmath>
#include <cstdio>
#include <cstdlib>

#include "constants.h"
#include "context.h"

int calc_shake_constraints(coord_t* coords, coord_t* xcoords) {
    auto& ctx = Context::instance();
    int total_iterations = 0;
    int shake = 0;

    for (int mol = 0; mol < ctx.n_molecules; mol++) {
        if (ctx.mol_n_shakes[mol] == 0) {
            continue;
        }
        int n_iterations = 0;

        bool converged = false;
        do {
            for (int i = 0; i < ctx.mol_n_shakes[mol]; i++) {
                ctx.shake_bonds[shake + i].ready = false;
            }

            for (int i = 0; i < ctx.mol_n_shakes[mol]; i++) {
                shake_bond_t& shake_bond = ctx.shake_bonds[shake + i];
                if (!shake_bond.ready) {
                    const int ai = shake_bond.ai - 1;
                    const int aj = shake_bond.aj - 1;
                    coord_t xij;
                    coord_t xxij;
                    double xij2, diff, corr, scp;

                    xij.x = coords[ai].x - coords[aj].x;
                    xij.y = coords[ai].y - coords[aj].y;
                    xij.z = coords[ai].z - coords[aj].z;
                    xij2 = std::pow(xij.x, 2) + std::pow(xij.y, 2) + std::pow(xij.z, 2);
                    diff = shake_bond.dist2 - xij2;
                    if (std::abs(diff) < shake_tol * shake_bond.dist2) {
                        shake_bond.ready = true;
                    }
                    xxij.x = xcoords[ai].x - xcoords[aj].x;
                    xxij.y = xcoords[ai].y - xcoords[aj].y;
                    xxij.z = xcoords[ai].z - xcoords[aj].z;
                    scp = xij.x * xxij.x + xij.y * xxij.y + xij.z * xxij.z;
                    corr = diff / (2.0 * scp * (ctx.winv[ai] + ctx.winv[aj]));

                    coords[ai].x += xxij.x * corr * ctx.winv[ai];
                    coords[ai].y += xxij.y * corr * ctx.winv[ai];
                    coords[ai].z += xxij.z * corr * ctx.winv[ai];
                    coords[aj].x -= xxij.x * corr * ctx.winv[aj];
                    coords[aj].y -= xxij.y * corr * ctx.winv[aj];
                    coords[aj].z -= xxij.z * corr * ctx.winv[aj];
                }
            }

            n_iterations++;

            converged = true;
            for (int i = 0; i < ctx.mol_n_shakes[mol]; i++) {
                if (!ctx.shake_bonds[shake + i].ready) {
                    converged = false;
                    break;
                }
            }
        } while (n_iterations < shake_max_iter && !converged);

        if (!converged) {
            for (int i = 0; i < ctx.mol_n_shakes[mol]; i++) {
                const int ai = ctx.shake_bonds[shake + i].ai - 1;
                const int aj = ctx.shake_bonds[shake + i].aj - 1;
                coord_t xxij;
                double xxij2;

                xxij.x = xcoords[ai].x - xcoords[aj].x;
                xxij.y = xcoords[ai].y - xcoords[aj].y;
                xxij.z = xcoords[ai].z - xcoords[aj].z;
                xxij2 = std::pow(xxij.x, 2) + std::pow(xxij.y, 2) + std::pow(xxij.z, 2);
                std::printf(">>> Shake failed, i = %d,j = %d, d = %f, d0 = %f", ai, aj, std::sqrt(xxij2),
                            ctx.shake_bonds[shake + i].dist2);
            }
            std::exit(EXIT_FAILURE);
        }

        shake += ctx.mol_n_shakes[mol];
        total_iterations += n_iterations;
    }

    // Set niter to the average number of iterations per molecule
    return ctx.n_molecules == 0 ? 0 : total_iterations / ctx.n_molecules;
}

void initial_shaking() {
    auto& ctx = Context::instance();

    for (int i = 0; i < ctx.n_atoms; i++) {
        ctx.xcoords[i].x = ctx.coords[i].x;
        ctx.xcoords[i].y = ctx.coords[i].y;
        ctx.xcoords[i].z = ctx.coords[i].z;
    }
    calc_shake_constraints(ctx.coords.data(), ctx.xcoords.data());
    for (int i = 0; i < ctx.n_atoms; i++) {
        ctx.xcoords[i].x = ctx.coords[i].x - ctx.dt * ctx.velocities[i].x;
        ctx.xcoords[i].y = ctx.coords[i].y - ctx.dt * ctx.velocities[i].y;
        ctx.xcoords[i].z = ctx.coords[i].z - ctx.dt * ctx.velocities[i].z;
    }
    calc_shake_constraints(ctx.xcoords.data(), ctx.coords.data());
    for (int i = 0; i < ctx.n_atoms; i++) {
        ctx.velocities[i].x = (ctx.coords[i].x - ctx.xcoords[i].x) / ctx.dt;
        ctx.velocities[i].y = (ctx.coords[i].y - ctx.xcoords[i].y) / ctx.dt;
        ctx.velocities[i].z = (ctx.coords[i].z - ctx.xcoords[i].z) / ctx.dt;
    }
}

void stop_cm_translation() {
    auto& ctx = Context::instance();
    auto &atypes = ctx.atypes->cpu_data_p;
    auto &catypes = ctx.catypes->cpu_data_p;
    double total_mass = 0;
    coord_t vcm = {};

    for (int ai = 0; ai < ctx.n_atoms; ai++) {
        const double rmass = catypes[atypes[ai].code - 1].m;
        total_mass += rmass;
        vcm.x += ctx.velocities[ai].x * rmass;
        vcm.y += ctx.velocities[ai].y;
        vcm.z += ctx.velocities[ai].z;
    }

    vcm.x /= total_mass;
    vcm.y /= total_mass;
    vcm.z /= total_mass;

    for (int ai = 0; ai < ctx.n_atoms; ai++) {
        ctx.velocities[ai].x -= vcm.x;
        ctx.velocities[ai].y -= vcm.y;
        ctx.velocities[ai].z -= vcm.z;
    }
}
