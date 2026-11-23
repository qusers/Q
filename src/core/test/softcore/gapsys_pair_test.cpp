// Unit tests for the Gapsys soft-core pair functions.
//
// Gapsys, Seeliger, de Groot (JCTC 2012): inside the soft-core radius rsc the
// pair potential is replaced by the second-order Taylor expansion of the hard
// potential around rsc,
//
//   V_soft(r) = V(rsc) + V'(rsc) (r - rsc) + 1/2 V''(rsc) (r - rsc)^2,   r < rsc
//
// with
//
//   rsc_lj = alpha_lj * (26/7 * sigma^6 * (1 - lambda))^(1/6)
//   rsc_q  = alpha_q * (1 + sigma_q * |qi * qj|) * (1 - lambda)^(1/6)
//
// The tests build that reference directly from these formulas and compare it
// with calc_vdw_gapsys / calc_electrostatic_gapsys (the kernels' polynomial
// form), apply_gapsys (non-bonded pairs) and calc_special_pair (1-4 pairs).
// Forces are checked against finite differences of the energy.
//
// Host only: no GPU is needed. The CPU and GPU kernels call these same functions.

#include <algorithm>
#include <cmath>
#include <cstdio>
#include <vector>

#include "common/include/bonded_force.h"
#include "common/include/pair_potential.h"

namespace {

int n_checks = 0;
int n_failures = 0;

// real_t is float in SPFP builds; the real_t functions are then only good to single precision.
const bool kDouble = sizeof(real_t) == sizeof(double);
const double kRealTol = kDouble ? 1.0e-11 : 2.0e-4;

void check_close(const char* what, double got, double want, double tol, double scale = 0.0) {
    n_checks++;
    const double s = std::max({1.0e-30, scale, std::fabs(want)});
    if (std::isfinite(got) && std::fabs(got - want) <= tol * s) return;
    n_failures++;
    std::printf("FAIL %s: got %.15g, want %.15g (tol %.1e, scale %.3g)\n", what, got, want, tol, s);
}

void check_true(const char* what, bool ok) {
    n_checks++;
    if (ok) return;
    n_failures++;
    std::printf("FAIL %s\n", what);
}

// ---- reference: the definition above, in double ----

struct Value {
    double v;   // energy
    double dv;  // dV/dr
};

// Second-order Taylor expansion around r0 of a function with value f, slope df and curvature d2f at r0.
Value taylor(double f, double df, double d2f, double r0, double r) {
    const double x = r - r0;
    return {f + df * x + 0.5 * d2f * x * x, df + d2f * x};
}

Value hard_lj(double A, double B, double r) {
    return {A / std::pow(r, 12) - B / std::pow(r, 6), -12.0 * A / std::pow(r, 13) + 6.0 * B / std::pow(r, 7)};
}

Value soft_lj(double A, double B, double r, double rsc) {
    if (r >= rsc) return hard_lj(A, B, r);
    const Value at = hard_lj(A, B, rsc);
    const double curvature = 156.0 * A / std::pow(rsc, 14) - 42.0 * B / std::pow(rsc, 8);
    return taylor(at.v, at.dv, curvature, rsc, r);
}

Value hard_coulomb(double qq, double r) {
    return {qq / r, -qq / (r * r)};
}

Value soft_coulomb(double qq, double r, double rsc) {
    if (r >= rsc) return hard_coulomb(qq, r);
    const Value at = hard_coulomb(qq, rsc);
    return taylor(at.v, at.dv, 2.0 * qq / std::pow(rsc, 3), rsc, r);
}

double rsc_lj(double alpha_lj, double sigma, double lambda) {
    return alpha_lj * std::pow(26.0 / 7.0 * std::pow(sigma, 6) * (1.0 - lambda), 1.0 / 6.0);
}

double rsc_q(double alpha_q, double sigma_q, double qi, double qj, double lambda) {
    return alpha_q * (1.0 + sigma_q * std::fabs(qi * qj)) * std::pow(1.0 - lambda, 1.0 / 6.0);
}

// ---- test data ----

constexpr double kCoulomb = 332.0716;

// One atom in the parameter form of each combination rule.
struct Atom {
    double a;  // arithmetic: R*; geometric: sqrt(C12)
    double b;  // arithmetic: sqrt(eps); geometric: sqrt(C6)
    double q;
};

struct PairCase {
    const char* name;
    int rule;
    Atom i, j;
};

const std::vector<PairCase> kPairs = {
    {"arithmetic C-O", VDW_ARITHMETIC, {1.908, std::sqrt(0.0860), 0.55}, {1.6612, std::sqrt(0.2100), -0.57}},
    {"arithmetic H-N", VDW_ARITHMETIC, {0.600, std::sqrt(0.0157), 0.27}, {1.824, std::sqrt(0.1700), -0.42}},
    {"arithmetic neutral", VDW_ARITHMETIC, {1.908, std::sqrt(0.1094), 0.0}, {1.908, std::sqrt(0.1094), 0.3}},
    {"geometric C-O", VDW_GEOMETRIC, {944.52, 22.03, 0.50}, {616.00, 23.25, -0.50}},
    {"geometric H-C", VDW_GEOMETRIC, {84.57, 5.42, 0.06}, {944.52, 22.03, -0.18}},
};

const std::vector<double> kLambdas = {0.0, 0.25, 0.5, 0.9, 0.999};

softcore_pair_params_t pair_params(double alpha_lj, double alpha_q, double sigma_q) {
    softcore_pair_params_t p;
    p.method = SC_GAPSYS;
    // Same expression as NonbondedForce::build_softcore_table.
    p.lj_scale = static_cast<real_t>(alpha_lj * std::pow(26.0 / 7.0, 1.0 / 6.0));
    p.alpha_q = static_cast<real_t>(alpha_q);
    p.sigma_q = static_cast<real_t>(sigma_q);
    return p;
}

// ---- tests ----

// The polynomial form equals the Taylor expansion, in double, over 0 <= r < rsc.
void test_polynomial_is_taylor_expansion() {
    const double A = 9.0e5, B = 6.0e2, qq = -0.3 * kCoulomb;
    for (double rsc : {0.4, 1.5, 3.1, 5.0}) {
        const Value lj_at = hard_lj(A, B, rsc);
        const Value q_at = hard_coulomb(qq, rsc);
        for (int k = 0; k <= 20; k++) {
            const double r = rsc * k / 20.0;

            const auto lj = calc_vdw_gapsys(Real2<double>{A, B}, r, rsc);
            const Value lj_ref = soft_lj(A, B, std::min(r, std::nextafter(rsc, 0.0)), rsc);
            // The expansion has terms of ~100 * V(rsc) that cancel; scale the tolerance by V(rsc).
            check_close("gapsys LJ energy", lj.x, lj_ref.v, 1e-11, 100.0 * (A / std::pow(rsc, 12) + B / std::pow(rsc, 6)));
            check_close("gapsys LJ dV/dr", lj.y, lj_ref.dv, 1e-11, 200.0 * (A / std::pow(rsc, 13) + B / std::pow(rsc, 7)));

            const auto el = calc_electrostatic_gapsys(-0.3, kCoulomb, r, rsc);
            const Value el_ref = soft_coulomb(qq, std::min(r, std::nextafter(rsc, 0.0)), rsc);
            check_close("gapsys Coulomb energy", el.x, el_ref.v, 1e-12, 3.0 * std::fabs(q_at.v));
            check_close("gapsys Coulomb dV/dr", el.y, el_ref.dv, 1e-12, 3.0 * std::fabs(q_at.dv));
        }

        // Value and slope are continuous at rsc.
        const auto lj = calc_vdw_gapsys(Real2<double>{A, B}, rsc, rsc);
        check_close("LJ continuous at rsc", lj.x, lj_at.v, 1e-11, 100.0 * (A / std::pow(rsc, 12) + B / std::pow(rsc, 6)));
        check_close("LJ slope continuous at rsc", lj.y, lj_at.dv, 1e-11, 200.0 * (A / std::pow(rsc, 13) + B / std::pow(rsc, 7)));
        const auto el = calc_electrostatic_gapsys(-0.3, kCoulomb, rsc, rsc);
        check_close("Coulomb continuous at rsc", el.x, q_at.v, 1e-12);
        check_close("Coulomb slope continuous at rsc", el.y, q_at.dv, 1e-12);
    }
}

// The energy at r = 0 is finite and, for a repulsive core, lower than the wall it replaces.
void test_finite_at_overlap() {
    const double A = 9.0e5, B = 6.0e2;
    const double rsc = 2.5;
    const auto lj = calc_vdw_gapsys(Real2<double>{A, B}, 0.0, rsc);
    const auto el = calc_electrostatic_gapsys(0.25, kCoulomb, 0.0, rsc);
    check_true("LJ finite at r = 0", std::isfinite(lj.x) && std::isfinite(lj.y));
    check_true("Coulomb finite at r = 0", std::isfinite(el.x) && std::isfinite(el.y));
    check_close("LJ at r = 0", lj.x, 91.0 * A / std::pow(rsc, 12) - 28.0 * B / std::pow(rsc, 6), 1e-12);
    check_close("Coulomb at r = 0", el.x, 3.0 * 0.25 * kCoulomb / rsc, 1e-12);
}

// sigma is where the hard LJ potential of the pair crosses zero.
void test_sigma_is_lj_zero_crossing() {
    for (const auto& c : kPairs) {
        const real_t sigma = softcore_pair_sigma(c.rule, c.i.a, c.i.b, c.j.a, c.j.b);
        const real_t2 pair = combine_vdw(c.rule, c.i.a, c.i.b, c.j.a, c.j.b);
        const real_t2 lj = calc_vdw(pair, static_cast<real_t>(1.0) / sigma);
        check_close(c.name, lj.x, 0.0, 20.0 * kRealTol, pair.x / std::pow(static_cast<double>(sigma), 12));
        check_true("sigma is positive", sigma > 0);
    }

    // Geometric rule without a zero crossing: 3 A if the pair has any LJ, 0 if it has none.
    check_close("geometric, C6 = 0", softcore_pair_sigma(VDW_GEOMETRIC, 900.0, 0.0, 600.0, 23.0), 3.0, kRealTol);
    check_close("geometric, C12 = 0", softcore_pair_sigma(VDW_GEOMETRIC, 0.0, 22.0, 600.0, 23.0), 3.0, kRealTol);
    check_close("geometric, no LJ", softcore_pair_sigma(VDW_GEOMETRIC, 0.0, 0.0, 600.0, 23.0), 0.0, kRealTol, 1.0);
}

// What the non-bonded kernels compute for one pair: hard potential, then apply_gapsys when sc > 0.
struct PairResult {
    double e_coul, e_vdw, dv;
};

PairResult kernel_pair(const PairCase& c, const softcore_pair_params_t& p, real_t sc, double r) {
    const real_t dis = static_cast<real_t>(r);
    const real_t inv_dis = static_cast<real_t>(1.0) / dis;
    const real_t qij = static_cast<real_t>(c.i.q) * static_cast<real_t>(c.j.q);
    const real_t2 pair = combine_vdw(c.rule, c.i.a, c.i.b, c.j.a, c.j.b);

    real_t2 el = calc_electrostatic(qij, kCoulomb, inv_dis);
    real_t2 vdw = calc_vdw(pair, inv_dis);
    if (sc > 0) {
        const real_t sigma = softcore_pair_sigma(c.rule, c.i.a, c.i.b, c.j.a, c.j.b);
        apply_gapsys(p, sc, sigma, qij, kCoulomb, pair, dis, el, vdw);
    }
    return {el.x, vdw.x, static_cast<double>(el.y) + vdw.y};
}

// apply_gapsys against the reference, for both rules, several lambdas, and r across both radii.
void test_nonbonded_pair_matches_reference() {
    const double alpha_lj = 0.85, alpha_q = 3.0, sigma_q = 1.0;
    const softcore_pair_params_t p = pair_params(alpha_lj, alpha_q, sigma_q);

    for (const auto& c : kPairs) {
        const real_t2 pair = combine_vdw(c.rule, c.i.a, c.i.b, c.j.a, c.j.b);
        const double A = pair.x, B = pair.y;
        const double qq = c.i.q * c.j.q * kCoulomb;
        // Independent of softcore_pair_sigma: the zero crossing of A / r^12 - B / r^6.
        const double sigma = std::pow(A / B, 1.0 / 6.0);

        for (double lambda : kLambdas) {
            const real_t sc = static_cast<real_t>(std::pow(1.0 - lambda, 1.0 / 6.0));
            const double r_lj = rsc_lj(alpha_lj, sigma, lambda);
            const double r_q = rsc_q(alpha_q, sigma_q, c.i.q, c.j.q, lambda);

            for (int k = 1; k <= 60; k++) {
                const double r = 0.1 * k;  // 0.1 .. 6.0 A
                // Skip points within rounding of a radius: either branch is right there.
                if (std::fabs(r - r_lj) < 1e-4 || std::fabs(r - r_q) < 1e-4) continue;

                const PairResult got = kernel_pair(c, p, sc, r);
                const Value lj = soft_lj(A, B, r, r_lj);
                const Value el = soft_coulomb(qq, r, r_q);

                const double r_eff = std::min(r, r_lj);
                const double lj_scale = 100.0 * (A / std::pow(r_eff, 12) + B / std::pow(r_eff, 6));
                check_close(c.name, got.e_vdw, lj.v, kRealTol, lj_scale);
                check_close(c.name, got.e_coul, el.v, kRealTol, 3.0 * std::fabs(qq) / std::min(r, r_q) + 1e-30);
                check_close(c.name, got.dv, lj.dv + el.dv, kRealTol, 2.0 * lj_scale / r_eff + 3.0 * std::fabs(qq) / (r * r));
            }
        }
    }
}

// lambda = 1 gives sc = 0: a pair of an enabled atom is then the hard pair, bit for bit.
// A pair outside both radii is also untouched.
void test_hard_limits() {
    const softcore_pair_params_t p = pair_params(0.85, 3.0, 1.0);
    for (const auto& c : kPairs) {
        for (int k = 5; k <= 60; k++) {
            const double r = 0.1 * k;
            const PairResult hard = kernel_pair(c, p, 0, r);

            const real_t sc_end = static_cast<real_t>(std::pow(1.0 - 1.0, 1.0 / 6.0));
            const PairResult end = kernel_pair(c, p, sc_end, r);
            check_true("lambda = 1 is the hard pair", end.e_vdw == hard.e_vdw && end.e_coul == hard.e_coul && end.dv == hard.dv);

            if (r > 9.0 * 0.85 && r > 3.0 * 2.0) {  // beyond any radius of these pairs
                const PairResult far = kernel_pair(c, p, 1, r);
                check_true("outside rsc is the hard pair", far.e_vdw == hard.e_vdw && far.e_coul == hard.e_coul && far.dv == hard.dv);
            }
        }
    }

    // A pair with no LJ and no charge stays zero, also at r inside the radii.
    const PairCase empty{"empty", VDW_GEOMETRIC, {0.0, 0.0, 0.0}, {600.0, 23.0, -0.5}};
    const PairResult z = kernel_pair(empty, p, 1, 0.5);
    check_true("pair without LJ or charge is zero", z.e_vdw == 0 && z.e_coul == 0 && z.dv == 0);
}

// dV/dr of the non-bonded pair equals the finite difference of its energy, also across the radii.
void test_nonbonded_force_is_energy_derivative() {
    if (!kDouble) return;  // a finite difference in single precision checks nothing useful
    const softcore_pair_params_t p = pair_params(0.85, 3.0, 1.0);
    const double h = 1.0e-5;
    for (const auto& c : kPairs) {
        for (double lambda : kLambdas) {
            const real_t sc = static_cast<real_t>(std::pow(1.0 - lambda, 1.0 / 6.0));
            for (int k = 2; k <= 120; k++) {
                const double r = 0.05 * k;
                const PairResult lo = kernel_pair(c, p, sc, r - h);
                const PairResult hi = kernel_pair(c, p, sc, r + h);
                const PairResult at = kernel_pair(c, p, sc, r);
                const double fd = ((hi.e_coul + hi.e_vdw) - (lo.e_coul + lo.e_vdw)) / (2.0 * h);
                // The curvature is continuous but the third derivative jumps at rsc; 1e-6 covers that.
                check_close(c.name, at.dv, fd, 1.0e-6, std::fabs(at.dv) + 1.0);
            }
        }
    }
}

// The energy has no jump when the radius sweeps over the pair as lambda changes.
void test_continuous_in_lambda() {
    if (!kDouble) return;
    const softcore_pair_params_t p = pair_params(0.85, 3.0, 1.0);
    const PairCase& c = kPairs[0];
    const double r = 1.8;
    double previous = 0;
    for (int k = 0; k <= 10000; k++) {
        const double lambda = k / 10000.0;
        const real_t sc = static_cast<real_t>(std::pow(1.0 - lambda, 1.0 / 6.0));
        const PairResult got = kernel_pair(c, p, sc, r);
        const double e = got.e_coul + got.e_vdw;
        // The steepest part is lambda -> 1, where sc = (1 - lambda)^(1/6) has an infinite slope.
        if (k > 0) check_true("energy continuous in lambda", std::fabs(e - previous) < 0.05 * (std::fabs(e) + 1.0));
        previous = e;
    }
}

// ---- 1-4 pairs ----

special_pair_term_t lj14_term(double A14, double B14, double qi, double qj, double el14_scale, double sigma,
                              double alpha_lj, double alpha_q, double sigma_q, bool soft) {
    special_pair_term_t term{};
    term.ai = 0;
    term.aj = 1;
    term.state = 0;
    term.kind = SpecialPairKind::LJ14;
    term.coulomb_coeff = kCoulomb * el14_scale * qi * qj;
    term.p1 = A14;
    term.p2 = B14;
    if (soft) {
        // Same expressions as build_lj14_pairs.
        term.sc_lj = alpha_lj * std::pow(26.0 / 7.0, 1.0 / 6.0) * sigma;
        term.sc_q = alpha_q * (1.0 + sigma_q * std::fabs(qi * qj));
    }
    return term;
}

void test_special_pair() {
    const double A = 4.5e5, B = 3.0e2, qi = 0.4, qj = -0.35, el14 = 0.8333;
    const double sigma = std::pow(A / B, 1.0 / 6.0);
    const double alpha_lj = 0.85, alpha_q = 3.0, sigma_q = 1.0;
    const special_pair_term_t soft = lj14_term(A, B, qi, qj, el14, sigma, alpha_lj, alpha_q, sigma_q, true);
    const special_pair_term_t hard = lj14_term(A, B, qi, qj, el14, sigma, alpha_lj, alpha_q, sigma_q, false);
    const double qq = kCoulomb * el14 * qi * qj;
    const double h = 1.0e-5;

    for (double lambda : kLambdas) {
        const double sc = std::pow(1.0 - lambda, 1.0 / 6.0);
        const double r_lj = rsc_lj(alpha_lj, sigma, lambda);
        // The 1-4 scale applies to the energy, not to the radius.
        const double r_q = rsc_q(alpha_q, sigma_q, qi, qj, lambda);

        for (int k = 2; k <= 120; k++) {
            const double r = 0.05 * k;
            if (std::fabs(r - r_lj) < 1e-4 || std::fabs(r - r_q) < 1e-4) continue;

            const auto got = calc_special_pair(soft, r, sc);
            const Value lj = soft_lj(A, B, r, r_lj);
            const Value el = soft_coulomb(qq, r, r_q);
            const double r_eff = std::min(r, r_lj);
            const double lj_scale = 100.0 * (A / std::pow(r_eff, 12) + B / std::pow(r_eff, 6));
            check_close("1-4 vdw", got.vdw, lj.v, 1e-11, lj_scale);
            check_close("1-4 coulomb", got.coul, el.v, 1e-11, 3.0 * std::fabs(qq) / std::min(r, r_q));
            check_close("1-4 dV/dr", got.dv, lj.dv + el.dv, 1e-11, 2.0 * lj_scale / r_eff + 3.0 * std::fabs(qq) / (r * r));

            const auto lo = calc_special_pair(soft, r - h, sc);
            const auto hi = calc_special_pair(soft, r + h, sc);
            const double fd = ((hi.coul + hi.vdw) - (lo.coul + lo.vdw)) / (2.0 * h);
            check_close("1-4 force is energy derivative", got.dv, fd, 1.0e-6, std::fabs(got.dv) + 1.0);
        }
    }

    for (int k = 10; k <= 120; k++) {
        const double r = 0.05 * k;
        const auto reference = calc_special_pair(hard, r);
        const Value lj = hard_lj(A, B, r);
        const Value el = hard_coulomb(qq, r);
        check_close("hard 1-4 vdw", reference.vdw, lj.v, 1e-12);
        check_close("hard 1-4 coulomb", reference.coul, el.v, 1e-12);
        check_close("hard 1-4 dV/dr", reference.dv, lj.dv + el.dv, 1e-12, std::fabs(lj.dv) + std::fabs(el.dv));

        // sc = 0 (lambda = 1, or a term that is not soft) is the hard pair, bit for bit.
        const auto end = calc_special_pair(soft, r, 0.0);
        check_true("1-4, lambda = 1 is the hard pair", end.vdw == reference.vdw && end.coul == reference.coul && end.dv == reference.dv);
        const auto unflagged = calc_special_pair(hard, r, 1.0);
        check_true("1-4 term without soft-core ignores sc", unflagged.vdw == reference.vdw && unflagged.coul == reference.coul && unflagged.dv == reference.dv);
    }

    // The exponential soft pair is a different kind; sc must not touch it.
    special_pair_term_t repulsion{};
    repulsion.kind = SpecialPairKind::SoftRepulsion;
    repulsion.coulomb_coeff = qq;
    repulsion.p1 = 40.0;
    repulsion.p2 = 2.5;
    const auto a = calc_special_pair(repulsion, 2.0);
    const auto b = calc_special_pair(repulsion, 2.0, 1.0);
    check_true("soft-repulsion pair ignores sc", a.vdw == b.vdw && a.coul == b.coul && a.dv == b.dv);
}

}  // namespace

int main() {
    test_polynomial_is_taylor_expansion();
    test_finite_at_overlap();
    test_sigma_is_lj_zero_crossing();
    test_nonbonded_pair_matches_reference();
    test_hard_limits();
    test_nonbonded_force_is_energy_derivative();
    test_continuous_in_lambda();
    test_special_pair();

    std::printf("%s precision: %d checks, %d failed\n", kDouble ? "double" : "single", n_checks, n_failures);
    return n_failures == 0 ? 0 : 1;
}
