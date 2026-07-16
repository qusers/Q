#pragma once
#include "context.h"
#include "temperature.h"

// Solvent-boundary model for the outer water shell around topo.solvent_center.
// Bundles the two terms that only act when solvent waters are present:
//   - radix: radial half-harmonic + Morse wall (needs R_TFREE for `shift`)
//   - polx:  boundary-water polarization  (only when md.polarisation)
struct WaterBoundaryData {
    // radix
    // U = Dwmz·(e^(2·awmz·db) − 2·e^(awmz·db))
    double Dwmz = 0;
    double awmz = 0;

    // polx
    int n_shells = 0;
    int n_max_inshell = 0;
    std::unique_ptr<HostDeviceBuffer<shell_t>> wshells;
    std::unique_ptr<HostDeviceBuffer<double>> theta;
    std::unique_ptr<HostDeviceBuffer<int>> list_sh;
};

class WaterBoundaryForce {
   public:
    virtual ~WaterBoundaryForce() = default;

    // Captures the temperature collaborator (radix reads R_TFREE for `shift`)
    // and runs backend setup (polx buffer allocation).
    void init(Context& ctx, const Temperature& temperature, const ParseResult& parsed);

    // One evaluation: radix always, polx when md.polarisation is on.
    // `iteration` drives polx's periodic theta_corr update (itdis_update).
    virtual void calc(Context& ctx, int iteration) = 0;

    virtual void cleanup() {}

    const WaterBoundaryData& data() const { return data_; }
    virtual void sync_for_output(Context& ctx);

   protected:
    WaterBoundaryData data_;
    const Temperature* temperature_ = nullptr;  // not owned; source of R_TFREE

    virtual void init_backend(Context& ctx) {}

   private:
    double build_water_sphere(Context& ctx);
    void build_wshells(Context& ctx, double crgQtot, const ParseResult& parsed);
};