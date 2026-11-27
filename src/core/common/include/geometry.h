#pragma once
#include <math.h>

#include "cuda_runtime_utility.h"
#include "md_types.h"

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

inline uint32_t expand_morton_bits(uint32_t value) {
    value &= 0x000003ffu;
    value = (value | (value << 16)) & 0x030000ffu;

    value = (value | (value << 8)) & 0x0300f00fu;

    value = (value | (value << 4)) & 0x030c30c3u;

    value = (value | (value << 2)) & 0x09249249u;

    return value;
}

inline uint32_t morton_code(uint32_t x, uint32_t y, uint32_t z) {
    return expand_morton_bits(x) | (expand_morton_bits(y) << 1) | (expand_morton_bits(z) << 2);
}
inline uint32_t quantize_morton_coordinate(double value, double minimum, double maximum) {
    const double extent = maximum - minimum;

    if (!(extent > 0.0)) {
        return 0;
    }

    double normalized = (value - minimum) / extent;
    normalized = std::max(0.0, std::min(1.0, normalized));

    return static_cast<uint32_t>(normalized * 1023.0 + 0.5);
}

inline uint32_t get_morton_code(const coord_t& p, double x_min, double x_max, double y_min, double y_max, double z_min, double z_max) {
    auto x = quantize_morton_coordinate(p.x, x_min, x_max);
    auto y = quantize_morton_coordinate(p.y, y_min, y_max);
    auto z = quantize_morton_coordinate(p.z, z_min, z_max);
    return morton_code(x, y, z);
}