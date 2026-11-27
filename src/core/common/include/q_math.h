#include "cuda_runtime_utility.h"
namespace {
HD inline double deg2rad(double d) {
    return d * (M_PI / 180.0);
}

double gauss(double mean, double sd) {
    double v1, v2, nd10;

    v1 = ((rand()) + 1.) / ((RAND_MAX) + 1.);
    v2 = ((rand()) + 1.) / ((RAND_MAX) + 1.);
    nd10 = cos(2 * M_PI * v2) * sqrt(-2. * log(v1));

    return sd * nd10 + mean;
}

}  // namespace
