#include <math.h>
#include <stdio.h>

#include "constants.h"
#include "common/include/precision.h"

// Get a value from a gaussian distributed random variable with
// mean mean and standard deviation sd
real_t gauss(real_t mean, real_t sd) {
    real_t v1, v2, nd10;

    v1 = ( (real_t)(rand()) + 1. )/( (real_t)(RAND_MAX) + 1. );
    v2 = ( (real_t)(rand()) + 1. )/( (real_t)(RAND_MAX) + 1. );
    nd10 = cos(2 * M_PI * v2) * sqrt(-2. * log(v1));

    return sd * nd10 + mean;
}


real_t to_degrees(real_t radians) {
    return radians * (180.0 / static_cast<real_t>(q_fortran_pi));
}

real_t to_radians(real_t degrees) {
    return degrees * (static_cast<real_t>(q_fortran_pi) / 180.0);
}
