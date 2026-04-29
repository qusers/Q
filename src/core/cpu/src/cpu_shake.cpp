#include "cpu_shake.h"

#include <cmath>
#include <cstdio>
#include <cstdlib>

#include "constants.h"
#include "context.h"

int calc_shake_constraints(coord_t* coords, coord_t* xcoords) {
    auto& ctx = Context::instance();
    auto *winv = ctx.winv->cpu_data_p;
    auto *mol_n_shakes = ctx.mol_n_shakes->cpu_data_p;
    auto *shake_bonds = ctx.shake_bonds->cpu_data_p;
    int total_iterations = 0;
    int shake = 0;

    for (int mol = 0; mol < ctx.n_molecules; mol++) {
        if (mol_n_shakes[mol] == 0) {
            continue;
        }
        int n_iterations = 0;

        bool converged = false;
        do {
            for (int i = 0; i < mol_n_shakes[mol]; i++) {
                shake_bonds[shake + i].ready = false;
            }

            for (int i = 0; i < mol_n_shakes[mol]; i++) {
                shake_bond_t& shake_bond = shake_bonds[shake + i];
                if (!shake_bond.ready) {
                    const int ai = shake_bond.ai - 1;
                    const int aj = shake_bond.aj - 1;
                    coord_t xij;
                    coord_t xxij;
                    real_t xij2, diff, corr, scp;

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
                    corr = diff / (2.0 * scp * (winv[ai] + winv[aj]));

                    coords[ai].x += xxij.x * corr * winv[ai];
                    coords[ai].y += xxij.y * corr * winv[ai];
                    coords[ai].z += xxij.z * corr * winv[ai];
                    coords[aj].x -= xxij.x * corr * winv[aj];
                    coords[aj].y -= xxij.y * corr * winv[aj];
                    coords[aj].z -= xxij.z * corr * winv[aj];
                }
            }

            n_iterations++;

            converged = true;
            for (int i = 0; i < mol_n_shakes[mol]; i++) {
                if (!shake_bonds[shake + i].ready) {
                    converged = false;
                    break;
                }
            }
        } while (n_iterations < shake_max_iter && !converged);

        if (!converged) {
            for (int i = 0; i < mol_n_shakes[mol]; i++) {
                const int ai = shake_bonds[shake + i].ai - 1;
                const int aj = shake_bonds[shake + i].aj - 1;
                coord_t xxij;
                real_t xxij2;

                xxij.x = xcoords[ai].x - xcoords[aj].x;
                xxij.y = xcoords[ai].y - xcoords[aj].y;
                xxij.z = xcoords[ai].z - xcoords[aj].z;
                xxij2 = std::pow(xxij.x, 2) + std::pow(xxij.y, 2) + std::pow(xxij.z, 2);
                std::printf(">>> Shake failed, i = %d,j = %d, d = %f, d0 = %f", ai, aj, std::sqrt(xxij2),
                            shake_bonds[shake + i].dist2);
            }
            std::exit(EXIT_FAILURE);
        }

        shake += mol_n_shakes[mol];
        total_iterations += n_iterations;
    }

    // Set niter to the average number of iterations per molecule
    return ctx.n_molecules == 0 ? 0 : total_iterations / ctx.n_molecules;
}

void initial_shaking() {
    auto& ctx = Context::instance();
    auto &coords = ctx.coords->cpu_data_p;
    auto &velocities = ctx.velocities->cpu_data_p;
    auto *xcoords = ctx.xcoords->cpu_data_p;

    for (int i = 0; i < ctx.n_atoms; i++) {
        xcoords[i].x = coords[i].x;
        xcoords[i].y = coords[i].y;
        xcoords[i].z = coords[i].z;
    }
    calc_shake_constraints(coords, xcoords);
    for (int i = 0; i < ctx.n_atoms; i++) {
        xcoords[i].x = coords[i].x - ctx.dt * velocities[i].x;
        xcoords[i].y = coords[i].y - ctx.dt * velocities[i].y;
        xcoords[i].z = coords[i].z - ctx.dt * velocities[i].z;
    }
    calc_shake_constraints(xcoords, coords);
    for (int i = 0; i < ctx.n_atoms; i++) {
        velocities[i].x = (coords[i].x - xcoords[i].x) / ctx.dt;
        velocities[i].y = (coords[i].y - xcoords[i].y) / ctx.dt;
        velocities[i].z = (coords[i].z - xcoords[i].z) / ctx.dt;
    }
}

void stop_cm_translation() {
    auto& ctx = Context::instance();
    auto &atypes = ctx.atypes->cpu_data_p;
    auto &catypes = ctx.catypes->cpu_data_p;
    auto &velocities = ctx.velocities->cpu_data_p;
    real_t total_mass = 0;
    coord_t vcm = {};

    for (int ai = 0; ai < ctx.n_atoms; ai++) {
        const real_t rmass = catypes[atypes[ai].code - 1].m;
        total_mass += rmass;
        vcm.x += velocities[ai].x * rmass;
        vcm.y += velocities[ai].y;
        vcm.z += velocities[ai].z;
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
