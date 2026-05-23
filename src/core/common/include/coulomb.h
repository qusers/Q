#pragma once

#include <math.h>

#include "precision.h"

template <typename Real>
QDYN_HOST_DEVICE inline Real fortran_scaled_charge(Real charge, Real coulomb_constant) {
    const Real rounded_charge = static_cast<Real>(static_cast<float>(charge));
    return static_cast<Real>(static_cast<float>(rounded_charge * sqrt(coulomb_constant)));
}

template <typename Real>
QDYN_HOST_DEVICE inline Real fortran_charge_product(Real qi, Real qj, Real coulomb_constant) {
    const Real qi_scaled = fortran_scaled_charge(qi, coulomb_constant);
    const Real qj_scaled = fortran_scaled_charge(qj, coulomb_constant);
    return static_cast<Real>(static_cast<float>(qi_scaled * qj_scaled));
}

template <typename Real>
QDYN_HOST_DEVICE inline Real fortran_q_charge_product(Real qi, Real qj, Real coulomb_constant) {
    const Real qi_scaled = fortran_scaled_charge(qi, coulomb_constant);
    const Real qj_scaled = fortran_scaled_charge(qj, coulomb_constant);
    return qi_scaled * qj_scaled;
}

template <typename Real>
QDYN_HOST_DEVICE inline Real fortran_coulomb_energy(
    Real qi,
    Real qj,
    Real rinv,
    Real coulomb_constant,
    Real scaling = static_cast<Real>(1.0)) {
    return fortran_charge_product(qi, qj, coulomb_constant) * rinv * scaling;
}

template <typename Real>
QDYN_HOST_DEVICE inline Real fortran_q_coulomb_energy(
    Real qi,
    Real qj,
    Real rinv,
    Real coulomb_constant,
    Real scaling = static_cast<Real>(1.0)) {
    return fortran_q_charge_product(qi, qj, coulomb_constant) * rinv * scaling;
}
