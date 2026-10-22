#pragma once

#include <cmath>

#include "constraint_force.h"
#include "context.h"
#include "host_device_buffer.h"
#include "precision.h"

struct bond_idx_t {
    int i, j;
};
struct angle_idx_t {
    int i, j, k;
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

struct special_pair_result_t {
    double coul;
    double vdw;
    double dv;
};

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
};
struct TorsionTerms {
    int n = 0;
    std::unique_ptr<HostDeviceBuffer<dihe_idx_t>> ids;
    std::unique_ptr<HostDeviceBuffer<torsion_param_t>> params;  // {k,n,d,paths} (real_t)
    std::unique_ptr<HostDeviceBuffer<int>> eslot;
    std::unique_ptr<HostDeviceBuffer<int>> state;  // -1 = ordinary topology term; >= 0 = Q-state index
};
struct ImproperTerms {  // improper2: harmonic, double
    int n = 0;
    std::unique_ptr<HostDeviceBuffer<dihe_idx_t>> ids;
    std::unique_ptr<HostDeviceBuffer<dparam2_t>> params;  // {k, phi0}   (double)
    std::unique_ptr<HostDeviceBuffer<int>> eslot;
    std::unique_ptr<HostDeviceBuffer<int>> state;  // -1 = ordinary topology term; >= 0 = Q-state index
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

HD inline dparam2_t calc_improper2(const double k, const double phi0, const double phi) {
    const double arg = 2.0 * phi - phi0;
    const double v = k * (1.0 + cos(arg));
    const double dv = k * -sin(arg) * 2.0;
    return {v, dv};
}

HD inline dparam2_t calc_q_improper(const double k, const double phi0, const double phi) {
    double delta = phi - phi0;
    const double two_pi = 2.0 * M_PI;
    delta -= two_pi * std::round(delta / two_pi);  // warps the angular difference into -pi to pi

    double v = 0.5 * k * delta * delta;
    double dv = k * delta;
    return {v, dv};
}

// Currently only SoftRepulsion reaches this function.
HD inline special_pair_result_t calc_special_pair(const special_pair_term_t& param, const double r) {
    const double inv_r = 1.0 / r;

    const double coul = param.coulomb_coeff * inv_r;
    const double vdw = param.p1 * exp(-param.p2 * r);

    const double dv = -coul * inv_r - param.p2 * vdw;

    return {coul, vdw, dv};
}