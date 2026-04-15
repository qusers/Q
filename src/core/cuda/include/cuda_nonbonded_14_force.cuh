#pragma once

#include <utility>

#include "common/include/nonbonded_14_mode.h"

void init_nonbonded_14_force_kernel_data();

void calc_nonbonded_14_forces_host();

void cleanup_nonbonded_14_force();
