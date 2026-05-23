#pragma once

#include "constants.h"
#include "precision.h"

QDYN_HOST_DEVICE inline real_t fortran_real(double value) {
    return static_cast<real_t>(static_cast<float>(value));
}

QDYN_HOST_DEVICE inline real_t fortran_radians_from_degrees(real_t degrees) {
    const double radians =
        (static_cast<double>(q_fortran_pi) / 180.0) * static_cast<double>(static_cast<float>(degrees));
    return static_cast<real_t>(static_cast<float>(radians));
}

QDYN_HOST_DEVICE inline real_t fortran_inverse_paths(real_t paths) {
    if (paths == static_cast<real_t>(0.0)) {
        return static_cast<real_t>(0.0);
    }
    return static_cast<real_t>(static_cast<float>(1.0f / static_cast<float>(paths)));
}
