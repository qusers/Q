#pragma once

#include "context.h"
#include "host_device_buffer.h"
#include "precision.h"

// Shared by CPU host build and CUDA; keep slot math in one place (like NB_HD).
#ifdef __CUDACC__
#define BF_HD __host__ __device__
#else
#define BF_HD
#endif

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

BF_HD inline int bonded_region_slots() { return 2; }
BF_HD inline int bonded_p_slot() { return 0; }
BF_HD inline int bonded_w_slot() { return 1; }

struct BondTerms {
    int n = 0;
    std::unique_ptr<HostDeviceBuffer<bond_idx_t>> ids;
    std::unique_ptr<HostDeviceBuffer<dparam2_t>> params;
    std::unique_ptr<HostDeviceBuffer<int>> eslot;
};
struct AngleTerms {
    int n = 0;
    std::unique_ptr<HostDeviceBuffer<angle_idx_t>> ids;
    std::unique_ptr<HostDeviceBuffer<dparam2_t>> params;  // {kth, th0} (double)
    std::unique_ptr<HostDeviceBuffer<int>> eslot;
};
struct TorsionTerms {
    int n = 0;
    std::unique_ptr<HostDeviceBuffer<dihe_idx_t>> ids;
    std::unique_ptr<HostDeviceBuffer<torsion_param_t>> params;  // {k,n,d,paths} (real_t)
    std::unique_ptr<HostDeviceBuffer<int>> eslot;
};
struct ImproperTerms {  // improper2: harmonic, double
    int n = 0;
    std::unique_ptr<HostDeviceBuffer<dihe_idx_t>> ids;
    std::unique_ptr<HostDeviceBuffer<dparam2_t>> params;  // {k, phi0}   (double)
    std::unique_ptr<HostDeviceBuffer<int>> eslot;
};

struct BondedData {
    BondTerms bond;
    AngleTerms angle;
    TorsionTerms torsion;
    ImproperTerms improper;

    bool enabled() const {
        return bond.n || angle.n || torsion.n || improper.n;
    }
};

class BondedForce {
   public:
    virtual ~BondedForce() = default;
    void init(Context& ctx);  // build_* then init_backend
    virtual void calc(Context& ctx) = 0;
    virtual void cleanup() {}
    bool enabled() const { return data_.enabled(); }
    BondedData& data() { return data_; }

   protected:
    BondedData data_;
    virtual void init_backend(Context& ctx) {}

   private:
    void build_bonds(Context& ctx);
    void build_angles(Context& ctx);
    void build_torsions(Context& ctx);
    void build_impropers(Context& ctx);
};

BF_HD inline dparam2_t calc_bond(const double k, const double r, const double r_eq) {
    // dv = dU / dr -> energy per length
    const double dr = r - r_eq;
    double v = 0.5 * k * dr * dr;
    double dv = k * dr;
    return {v, dv};
}

BF_HD inline dparam2_t calc_angle(const double k, const double th, const double th_eq) {
    // dv = dU / dth. -> energy per angle
    const double dth = th - th_eq;
    double v = 0.5 * k * dth * dth;
    double dv = k * dth;
    return {v, dv};
}

BF_HD inline real_t2 calc_torsion(const real_t k, const int n, const real_t phi, const real_t gamma, real_t paths) {
    // dv = dU / dphi -> energy per angle
    const real_t arg = n * phi - gamma;
    const real_t v = k * (1 + cos(arg)) * paths;
    const real_t dv = -k * n * sin(arg) * paths;
    return {v, dv};
}

BF_HD inline dparam2_t calc_improper2(const double k, const double phi0, const double phi) {
    const double arg = 2.0 * phi - phi0;
    const double v = k * (1.0 + cos(arg));
    const double dv = k * -sin(arg) * 2.0;
    return {v, dv};
}