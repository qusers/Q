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

/* =============================================
 * == CALCUTED IN THE INTEGRATION
 * =============================================
 */

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
