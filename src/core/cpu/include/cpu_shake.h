#pragma once

#include "md_types.h"

class Context;

int calc_shake_constraints(Context& ctx, coord_t* coords, coord_t* xcoords);
void initial_shaking(Context& ctx);
void stop_cm_translation(Context& ctx);
