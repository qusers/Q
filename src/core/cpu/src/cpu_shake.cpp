#include "cpu_shake.h"

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
        const double dt = static_cast<double>(ctx.dt);
        xcoords[i].x = static_cast<real_t>(static_cast<double>(coords[i].x) - dt * static_cast<double>(velocities[i].x));
        xcoords[i].y = static_cast<real_t>(static_cast<double>(coords[i].y) - dt * static_cast<double>(velocities[i].y));
        xcoords[i].z = static_cast<real_t>(static_cast<double>(coords[i].z) - dt * static_cast<double>(velocities[i].z));
    }

    apply_to(ctx, xcoords, coords);

    for (int i = 0; i < ctx.n_atoms; i++) {
        const double dt = static_cast<double>(ctx.dt);
        velocities[i].x = static_cast<real_t>((static_cast<double>(coords[i].x) - static_cast<double>(xcoords[i].x)) / dt);
        velocities[i].y = static_cast<real_t>((static_cast<double>(coords[i].y) - static_cast<double>(xcoords[i].y)) / dt);
        velocities[i].z = static_cast<real_t>((static_cast<double>(coords[i].z) - static_cast<double>(xcoords[i].z)) / dt);
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

                const double xij_x = static_cast<double>(coords[ai].x) - static_cast<double>(coords[aj].x);
                const double xij_y = static_cast<double>(coords[ai].y) - static_cast<double>(coords[aj].y);
                const double xij_z = static_cast<double>(coords[ai].z) - static_cast<double>(coords[aj].z);

                const double current_dist2 = xij_x * xij_x + xij_y * xij_y + xij_z * xij_z;
                const double diff = shake_bonds[shake + i].dist2 - current_dist2;

                if (std::abs(diff) < shake_tol * shake_bonds[shake + i].dist2) {
                    ready[shake + i] = true;
                } else {
                    converged = false;
                }

                const double xxij_x = static_cast<double>(xcoords[ai].x) - static_cast<double>(xcoords[aj].x);
                const double xxij_y = static_cast<double>(xcoords[ai].y) - static_cast<double>(xcoords[aj].y);
                const double xxij_z = static_cast<double>(xcoords[ai].z) - static_cast<double>(xcoords[aj].z);
                const double scp = xij_x * xxij_x + xij_y * xxij_y + xij_z * xxij_z;
                const double winv_ai = static_cast<double>(winv[ai]);
                const double winv_aj = static_cast<double>(winv[aj]);
                const double corr = diff / (2.0 * scp * (winv_ai + winv_aj));

                coords[ai].x = static_cast<real_t>(static_cast<double>(coords[ai].x) + xxij_x * corr * winv_ai);
                coords[ai].y = static_cast<real_t>(static_cast<double>(coords[ai].y) + xxij_y * corr * winv_ai);
                coords[ai].z = static_cast<real_t>(static_cast<double>(coords[ai].z) + xxij_z * corr * winv_ai);

                coords[aj].x = static_cast<real_t>(static_cast<double>(coords[aj].x) - xxij_x * corr * winv_aj);
                coords[aj].y = static_cast<real_t>(static_cast<double>(coords[aj].y) - xxij_y * corr * winv_aj);
                coords[aj].z = static_cast<real_t>(static_cast<double>(coords[aj].z) - xxij_z * corr * winv_aj);
            }

            try_times++;
        } while (try_times < shake_max_iter && !converged);

        if (!converged) {
            for (int i = 0; i < current_shake_num; i++) {
                const int ai = shake_bonds[shake + i].ai - 1;
                const int aj = shake_bonds[shake + i].aj - 1;

                const double xxij_x = static_cast<double>(xcoords[ai].x) - static_cast<double>(xcoords[aj].x);
                const double xxij_y = static_cast<double>(xcoords[ai].y) - static_cast<double>(xcoords[aj].y);
                const double xxij_z = static_cast<double>(xcoords[ai].z) - static_cast<double>(xcoords[aj].z);
                const double dist2 = xxij_x * xxij_x + xxij_y * xxij_y + xxij_z * xxij_z;

                printf(">>> Shake failed, i = %d, j = %d, d = %f, d0 = %f", ai, aj, sqrt(dist2), sqrt(shake_bonds[shake + i].dist2));
            }
            std::exit(EXIT_FAILURE);
        }

        shake += current_shake_num;
    }
}
