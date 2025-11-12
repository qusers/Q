#ifndef __UTILS_H__
#define __UTILS_H__

/* =============================================
 * == MATH STUFF
 * =============================================
 */

float gauss(float mean, float sd);
float to_degrees(float radians);
float to_radians(float degrees);

/* =============================================
 * == DEVICE
 * =============================================
 */

void check_cudaMalloc(void** devPtr, size_t size);

#endif /* __UTILS_H__ */