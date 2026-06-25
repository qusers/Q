#pragma once
#include <math.h>

#include "md_types.h"

#ifdef __CUDACC__
#define GEO_HD __host__ __device__
#else
#define GEO_HD
#endif

template <class T>
GEO_HD inline Real3<T> operator-(const Real3<T>& a, const Real3<T>& b) {
    return {a.x - b.x, a.y - b.y, a.z - b.z};
}
template <class T>
GEO_HD inline Real3<T> operator+(const Real3<T>& a, const Real3<T>& b) {
    return {a.x + b.x, a.y + b.y, a.z + b.z};
}
template <class T, class S>
GEO_HD inline Real3<T> operator*(const Real3<T>& a, S s) {
    const T t = static_cast<T>(s);
    return {a.x * t, a.y * t, a.z * t};
}
template <class T, class S>
GEO_HD inline Real3<T> operator*(S s, const Real3<T>& a) { return a * s; }
template <class T, class S>
GEO_HD inline Real3<T> operator/(const Real3<T>& a, S s) {
    const T t = static_cast<T>(s);
    return {a.x / t, a.y / t, a.z / t};
}

template <class T>
GEO_HD inline T dot(const Real3<T>& a, const Real3<T>& b) {
    return a.x * b.x + a.y * b.y + a.z * b.z;
}
template <class T>
GEO_HD inline Real3<T> cross(const Real3<T>& a, const Real3<T>& b) {
    return {a.y * b.z - a.z * b.y,
            a.z * b.x - a.x * b.z,
            a.x * b.y - a.y * b.x};
}
template <class T>
GEO_HD inline T norm2(const Real3<T>& a) { return dot(a, a); }  // squared length

template <class T>
GEO_HD inline T norm(const Real3<T>& a) { return sqrt(norm2(a)); }  // length

template <class T>
GEO_HD inline Real3<T> get_perpendicular_vector(const Real3<T> a, const Real3<T>& b) {
    // Get the vector that is perpendicular to `a` inside the plane formed by vector a, b
    Real3<T> a_unit = a / norm(a);
    Real3<T> b_unit = b / norm(b);
    double cos_th = dot(a_unit, b_unit);
    return b_unit - a_unit * cos_th;
}
