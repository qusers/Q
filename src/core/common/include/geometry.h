#pragma once
#include <math.h>

#include "md_types.h"
#include "cuda_runtime_utility.h"


template <class T>
HD inline Real3<T> operator-(const Real3<T>& a, const Real3<T>& b) {
    return {a.x - b.x, a.y - b.y, a.z - b.z};
}
template <class T>
HD inline Real3<T> operator+(const Real3<T>& a, const Real3<T>& b) {
    return {a.x + b.x, a.y + b.y, a.z + b.z};
}
template <class T, class S>
HD inline Real3<T> operator*(const Real3<T>& a, S s) {
    const T t = static_cast<T>(s);
    return {a.x * t, a.y * t, a.z * t};
}
template <class T, class S>
HD inline Real3<T> operator*(S s, const Real3<T>& a) { return a * s; }
template <class T, class S>
HD inline Real3<T> operator/(const Real3<T>& a, S s) {
    const T t = static_cast<T>(s);
    return {a.x / t, a.y / t, a.z / t};
}

template <class T>
HD inline T dot(const Real3<T>& a, const Real3<T>& b) {
    return a.x * b.x + a.y * b.y + a.z * b.z;
}
template <class T>
HD inline Real3<T> cross(const Real3<T>& a, const Real3<T>& b) {
    return {a.y * b.z - a.z * b.y,
            a.z * b.x - a.x * b.z,
            a.x * b.y - a.y * b.x};
}
template <class T>
HD inline T norm2(const Real3<T>& a) { return dot(a, a); }  // squared length

template <class T>
HD inline T norm(const Real3<T>& a) { return sqrt(norm2(a)); }  // length

template <class T>
HD inline Real3<T> get_perpendicular_vector(const Real3<T> a, const Real3<T>& b) {
    // Get the vector that is perpendicular to `a` inside the plane formed by vector a, b
    Real3<T> a_unit = a / norm(a);
    Real3<T> b_unit = b / norm(b);
    T cos_th = dot(a_unit, b_unit);
    return b_unit - a_unit * cos_th;
}

template <class To, class From>
HD inline Real3<To> real3_cast(const Real3<From>& a) {
    return {
        static_cast<To>(a.x),
        static_cast<To>(a.y),
        static_cast<To>(a.z)};
}
