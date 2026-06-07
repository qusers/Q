#include "cpu_shake_v2.h"

#include <cmath>
#include <vector>

#include "constants.h"

void CpuShake::apply(Context& ctx) {
    auto* coords = ctx.coords->cpu_data_p;
    auto* xcoords = ctx.xcoords->cpu_data_p;
    apply_to(ctx, coords, xcoords);
}

void CpuShake::initial_shake(Context& ctx) {
    /*
     * Fresh-start setup has no previous coordinate frame yet.
     *
     * The leapfrog/SHAKE state uses:
     *   coords  = current coordinates
     *   xcoords = previous/reference coordinates
     *   velocity = (coords - xcoords) / dt
     *
     * Initial SHAKE should therefore:
     *   1. Copy coords into xcoords and SHAKE coords so the initial structure
     *      satisfies all distance constraints.
     *   2. Reconstruct a previous frame with xcoords = coords - dt * velocity.
     *   3. SHAKE that reconstructed previous frame as well, using coords as
     *      the reference frame.
     *   4. Recompute velocities from the corrected coords/xcoords pair so the
     *      initial velocities are consistent with the constrained geometry.
     */
    auto* coords = ctx.coords->cpu_data_p;
    auto* velocities = ctx.velocities->cpu_data_p;
    auto* xcoords = ctx.xcoords->cpu_data_p;

    for (int i = 0; i < ctx.n_atoms; i++) {
        xcoords[i] = coords[i];
    }

    apply_to(ctx, coords, xcoords);

    for (int i = 0; i < ctx.n_atoms; i++) {
        xcoords[i].x = coords[i].x - ctx.dt * velocities[i].x;
        xcoords[i].y = coords[i].y - ctx.dt * velocities[i].y;
        xcoords[i].z = coords[i].z - ctx.dt * velocities[i].z;
    }

    apply_to(ctx, xcoords, coords);

    for (int i = 0; i < ctx.n_atoms; i++) {
        velocities[i].x = (coords[i].x - xcoords[i].x) / ctx.dt;
        velocities[i].y = (coords[i].y - xcoords[i].y) / ctx.dt;
        velocities[i].z = (coords[i].z - xcoords[i].z) / ctx.dt;
    }
}

void CpuShake::apply_to(Context& ctx, coord_t* coords, coord_t* xcoords) {
    /*
     * This follows the original Fortran SHAKE loop: each molecule owns a
     * ready flag per constraint, and a constraint that reaches tolerance is
     * skipped for the rest of that molecule's iteration loop.
     *
     * That behavior is kept for parity with Q, but it is not obviously the
     * safest iterative scheme. If two constraints share atoms, a later
     * correction can move a previously-ready bond away from tolerance again.
     * A stricter implementation would recheck every constraint each iteration,
     * or at least clear ready after coupled corrections.
     */
    const auto* winv = ctx.winv->cpu_data_p;
    const auto* shake_bonds = data_.shake_bonds->cpu_data_p;
    const auto* mol_n_shakes = data_.mol_n_shakes->cpu_data_p;

    // Kept outside the iteration loop to match the Fortran ready semantics.
    std::vector<bool> ready(data_.n_constraints, false);
    int shake = 0;
    for (int mol = 0; mol < ctx.n_molecules; mol++) {
        if (mol_n_shakes[mol] == 0) continue;

        bool converged;
        int try_times = 0;

        int current_shake_num = mol_n_shakes[mol];
        do {
            converged = true;
            for (int i = 0; i < current_shake_num; i++) {
                if (ready[shake + i]) continue;

                const int ai = shake_bonds[shake + i].ai - 1;
                const int aj = shake_bonds[shake + i].aj - 1;

                coord_t xij = {coords[ai].x - coords[aj].x, coords[ai].y - coords[aj].y, coords[ai].z - coords[aj].z};

                real_t current_dist2 = xij.x * xij.x + xij.y * xij.y + xij.z * xij.z;
                real_t diff = shake_bonds[shake + i].dist2 - current_dist2;

                if (std::abs(diff) < shake_tol * shake_bonds[shake + i].dist2) {
                    ready[shake + i] = true;
                } else {
                    converged = false;
                }

                coord_t xxij = {xcoords[ai].x - xcoords[aj].x, xcoords[ai].y - xcoords[aj].y,
                                xcoords[ai].z - xcoords[aj].z};
                real_t scp = xij.x * xxij.x + xij.y * xxij.y + xij.z * xxij.z;
                real_t corr = diff / (2.0 * scp * (winv[ai] + winv[aj]));

                coords[ai].x += xxij.x * corr * winv[ai];
                coords[ai].y += xxij.y * corr * winv[ai];
                coords[ai].z += xxij.z * corr * winv[ai];

                coords[aj].x -= xxij.x * corr * winv[aj];
                coords[aj].y -= xxij.y * corr * winv[aj];
                coords[aj].z -= xxij.z * corr * winv[aj];
            }

            try_times++;

        } while (try_times < shake_max_iter && !converged);

        if (!converged) {
            for (int i = 0; i < current_shake_num; i++) {
                const int ai = shake_bonds[shake + i].ai - 1;
                const int aj = shake_bonds[shake + i].aj - 1;

                coord_t xxij = {xcoords[ai].x - xcoords[aj].x, xcoords[ai].y - xcoords[aj].y, xcoords[ai].z - xcoords[aj].z};
                real_t dist2 = xxij.x * xxij.x + xxij.y * xxij.y + xxij.z * xxij.z;

                printf(">>> Shake failed, i = %d, j = %d, d^2 = %f, d0^2 = %f", ai, aj, dist2, shake_bonds[shake + i].dist2);
            }
            std::exit(EXIT_FAILURE);
        }

        shake += current_shake_num;
    }
}
