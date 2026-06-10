#pragma once
#include "shake.h"

void init_leapfrog_kernel_data();
void calc_leapfrog_host(Shake& shake);
void cleanup_leapfrog();
