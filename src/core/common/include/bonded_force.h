#pragma once

#include <cmath>

#include "constraint_force.h"
#include "context.h"
#include "host_device_buffer.h"
#include "pair_potential.h"
#include "precision.h"

struct bond_idx_t {
    int i, j;
};
struct angle_idx_t {
    int i, j, k;
};

struct morse_coupling_t {
    int i = -1;
    int j = -1;
    double a = 0.0;
    double r0 = 0.0;
    bool inverted = false;
};

struct dihe_idx_t {
    int i, j, k, l;
};  // torsion + improper2

struct torsion_param_t {
    real_t k, n, d, paths;
};  // working precision (PMEFloat)

struct dparam2_t {
    double x, y;
};

template <typename T>
struct SpecialPairResult {
    T coul;
    T vdw;
    T dv;
};
using special_pair_result_t = SpecialPairResult<double>;

struct SpecialPairTerms {
    int n = 0;
    std::unique_ptr<HostDeviceBuffer<special_pair_term_t>> terms;
};

struct BondTerms {
    int n = 0;
    std::unique_ptr<HostDeviceBuffer<bond_idx_t>> ids;
    std::unique_ptr<HostDeviceBuffer<cbond_t>> params;
    std::unique_ptr<HostDeviceBuffer<int>> eslot;
    std::unique_ptr<HostDeviceBuffer<int>> state;  // -1 = ordinary topology term; >= 0 = Q-state index
};
struct AngleTerms {
    int n = 0;
    std::unique_ptr<HostDeviceBuffer<angle_idx_t>> ids;
    std::unique_ptr<HostDeviceBuffer<dparam2_t>> params;  // {kth, th0} (double)
    std::unique_ptr<HostDeviceBuffer<int>> eslot;
    std::unique_ptr<HostDeviceBuffer<int>> state;  // -1 = ordinary topology term; >= 0 = Q-state index

    std::unique_ptr<HostDeviceBuffer<morse_coupling_t>> couplings;
};
struct TorsionTerms {
    int n = 0;
    std::unique_ptr<HostDeviceBuffer<dihe_idx_t>> ids;
    std::unique_ptr<HostDeviceBuffer<torsion_param_t>> params;  // {k,n,d,paths} (real_t)
    std::unique_ptr<HostDeviceBuffer<int>> eslot;
    std::unique_ptr<HostDeviceBuffer<int>> state;  // -1 = ordinary topology term; >= 0 = Q-state index

    std::unique_ptr<HostDeviceBuffer<morse_coupling_t>> couplings;
};
struct ImproperTerms {  // improper2: harmonic, double
    int n = 0;
    std::unique_ptr<HostDeviceBuffer<dihe_idx_t>> ids;
    std::unique_ptr<HostDeviceBuffer<dparam2_t>> params;  // {k, phi0}   (double)
    std::unique_ptr<HostDeviceBuffer<int>> eslot;
    std::unique_ptr<HostDeviceBuffer<int>> state;  // -1 = ordinary topology term; >= 0 = Q-state index

    std::unique_ptr<HostDeviceBuffer<morse_coupling_t>> couplings;
};

struct BondedData {
    BondTerms bond;
    AngleTerms angle;
    TorsionTerms torsion;
    ImproperTerms improper;
    SpecialPairTerms special_pairs;

    bool enabled() const {
        return bond.n || angle.n || torsion.n || improper.n || special_pairs.n;
    }
};

class BondedForce {
   public:
    virtual ~BondedForce() = default;
    void init(Context& ctx, const ParseResult& parsed, const ConstraintData& constraint_data, const std::vector<special_pair_term_t>& special_pairs = {});  // build_* then init_backend
    virtual void calc(Context& ctx) = 0;
    virtual void cleanup() {}
    bool enabled() const { return data_.enabled(); }
    BondedData& data() { return data_; }

   protected:
    BondedData data_;
    virtual void init_backend(Context& ctx) {}

   private:
    void build_bonds(Context& ctx, const ParseResult& parsed, const ConstraintData& shake_data);
    void build_angles(Context& ctx, const ParseResult& parsed, const ConstraintData& shake_data);
    void build_torsions(Context& ctx, const ParseResult& parsed, const ConstraintData& shake_data);
    void build_impropers(Context& ctx, const ParseResult& parsed, const ConstraintData& shake_data);
    void build_special_pairs(Context& ctx, const std::vector<special_pair_term_t>& special_pairs);
};

// TorsionTerms::state of a topology torsion whose four atoms are Q atoms: weight 1 like any
// topology term, times the REST2 scale.
constexpr int REST_TORSION = -2;

HD inline double special_pair_rest(const special_pair_term_t& term, double rest, double rest_sqrt) {
    return term.n_q == 2 ? rest : term.n_q == 1 ? rest_sqrt
                                                : 1.0;
}

HD inline dparam2_t calc_morse_bond(const double D, const double a, const double r, const double r_eq) {
    const double dr = r - r_eq;
    const double f = exp(-a * dr);

    const double v = D * (f * f - 2.0 * f);

    const double dv = 2.0 * D * a * (f - f * f);
    return {v, dv};
}

HD inline dparam2_t calc_harmonic_bond(const double k, const double r, const double r_eq) {
    // dv = dU / dr -> energy per length
    const double dr = r - r_eq;
    double v = 0.5 * k * dr * dr;
    double dv = k * dr;
    return {v, dv};
}

HD inline dparam2_t calc_bond(const cbond_t& param, const double r) {
    if (!param.is_morse) {
        return calc_harmonic_bond(param.kb, r, param.b0);
    } else {
        return calc_morse_bond(param.D, param.a, r, param.b0);
    }
}

// Returns {gamma, dgamma/dr}
HD inline dparam2_t calc_morse_coupling(const double a, const double r, const double r0) {
    if (r <= r0) {
        return {1.0, 0.0};
    }

    const double f = exp(-a * (r - r0));
    const double gamma = 2.0 * f - f * f;
    const double dgamma_dr = -2.0 * a * (f - f * f);

    return {gamma, dgamma_dr};
}

HD inline dparam2_t calc_angle(const double k, const double th, const double th_eq) {
    // dv = dU / dth. -> energy per angle
    const double dth = th - th_eq;
    double v = 0.5 * k * dth * dth;
    double dv = k * dth;
    return {v, dv};
}

HD inline real_t2 calc_torsion(const real_t k, const real_t n, const real_t phi, const real_t gamma, real_t paths) {
    // dv = dU / dphi -> energy per angle
    const real_t arg = n * phi - gamma;
    const real_t v = k * (1 + cos(arg)) * paths;
    const real_t dv = -k * n * sin(arg) * paths;
    return {v, dv};
}

// T is the precision of the geometry: double, or real_t in the GPU kernel.
template <typename T>
HD inline Real2<T> calc_improper2(const T k, const T phi0, const T phi) {
    const T arg = static_cast<T>(2) * phi - phi0;
    const T v = k * (static_cast<T>(1) + cos(arg));
    const T dv = k * -sin(arg) * static_cast<T>(2);
    return {v, dv};
}

template <typename T>
HD inline Real2<T> calc_q_improper(const T k, const T phi0, const T phi) {
    T delta = phi - phi0;
    const T two_pi = static_cast<T>(2.0 * M_PI);
    delta -= two_pi * round(delta / two_pi);  // warps the angular difference into -pi to pi

    T v = static_cast<T>(0.5) * k * delta * delta;
    T dv = k * delta;
    return {v, dv};
}

// sc = (1 - lambda)^(1/6) for a soft 1-4 pair, 0 otherwise.
// T is the precision of the geometry: double, or real_t in the GPU kernel.
template <typename T>
HD inline SpecialPairResult<T> calc_special_pair(const special_pair_term_t& param, const T r, const T sc = 0) {
    const T inv_r = static_cast<T>(1) / r;

    T coul = static_cast<T>(param.coulomb_coeff) * inv_r;
    T dv_coul = -coul * inv_r;

    T vdw, dv_vdw;

    if (param.kind == SpecialPairKind::LJ14) {
        const T inv_r2 = inv_r * inv_r;
        const T inv_r6 = inv_r2 * inv_r2 * inv_r2;

        const T Va = static_cast<T>(param.p1) * inv_r6 * inv_r6;
        const T Vb = static_cast<T>(param.p2) * inv_r6;

        vdw = Va - Vb;
        dv_vdw = (static_cast<T>(-12) * Va + static_cast<T>(6) * Vb) * inv_r;

        if (sc > 0) {
            const T rsc_lj = static_cast<T>(param.sc_lj) * sc;
            if (r < rsc_lj) {
                const auto soft = calc_vdw_gapsys(Real2<T>{static_cast<T>(param.p1), static_cast<T>(param.p2)}, r, rsc_lj);
                vdw = soft.x;
                dv_vdw = soft.y;
            }

            const T rsc_q = static_cast<T>(param.sc_q) * sc;
            if (r < rsc_q) {
                // coulomb_coeff already holds the Coulomb constant and the 1-4 scale.
                const auto soft = calc_electrostatic_gapsys(static_cast<T>(param.coulomb_coeff), static_cast<T>(1), r, rsc_q);
                coul = soft.x;
                dv_coul = soft.y;
            }
        }
    } else {
        vdw = static_cast<T>(param.p1) * exp(-static_cast<T>(param.p2) * r);
        dv_vdw = -static_cast<T>(param.p2) * vdw;
    }
    const T dv = dv_coul + dv_vdw;

    return {coul, vdw, dv};
}