#include <math.h>
#include <stdio.h>

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
    return radians * (180.0 / M_PI);
}




real_t to_radians(real_t degrees) {
    const real_t PI = 4.0 * atanf(1.0);
    return degrees * (PI / 180.0);
}
