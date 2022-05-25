#include <math.h>
#include <stdio.h>

// Get a value from a gaussian distributed random variable with
// mean mean and standard deviation sd
double gauss(double mean, double sd) {
    double v1, v2, nd10;

    v1 = ( (double)(rand()) + 1. )/( (double)(RAND_MAX) + 1. );
    v2 = ( (double)(rand()) + 1. )/( (double)(RAND_MAX) + 1. );
    nd10 = cos(2 * M_PI * v2) * sqrt(-2. * log(v1));

    return sd * nd10 + mean;
}


double to_degrees(double radians) {
    return radians * (180.0 / M_PI);
}

double to_radians(double degrees) {
    return degrees * (M_PI / 180.0);
}

void check_cudaMalloc(void** devPtr, size_t size) {
    cudaError_t error = cudaMalloc((void**) devPtr, size);
    if (error != cudaSuccess) {
        printf(">>> FATAL: memory for matrix could not be allocated. Exiting...\n");
        exit(EXIT_FAILURE);
    }
}