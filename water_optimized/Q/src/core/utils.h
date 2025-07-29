#ifndef __UTILS_H__
#define __UTILS_H__

/* =============================================
 * == MATH STUFF
 * =============================================
 */

double gauss(double mean, double sd);
double to_degrees(double radians);
double to_radians(double degrees);

/* =============================================
 * == DEVICE
 * =============================================
 */

void check_cudaMalloc(void** devPtr, size_t size);

#endif /* __UTILS_H__ */