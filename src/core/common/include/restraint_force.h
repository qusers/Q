#pragma once
#include "context.h"
#include "host_device_buffer.h"
#include "precision.h"

struct PosRestraintTerms {
    int n = 0;
    std::unique_ptr<HostDeviceBuffer<int>> ids;      // atom index (0-based)
    std::unique_ptr<HostDeviceBuffer<double>> k;     // k_fix or k_pshell
    std::unique_ptr<HostDeviceBuffer<int>> eslot_a;  // E_RESTR_FIX  or -1
    std::unique_ptr<HostDeviceBuffer<int>> eslot_b;  // E_RESTR_SHELL or -1
};

// --- restrpos: explicit ref point, anisotropic k, FEP-aware (ipsi) ---
struct AtomRestraintTerms {
    int n = 0;
    std::unique_ptr<HostDeviceBuffer<restrpos_t>> recs;  // {a, ipsi, x(ref), k(kx,ky,kz)}
};

// --- restrseq: group range [ai,aj), 3 modes, ref = coords_init ---
struct SeqRestraintTerms {
    int n = 0;
    std::unique_ptr<HostDeviceBuffer<restrseq_t>> recs;  // {ai, aj, k, ih, to_center}
};

// --- restrdis: atom pair, flat-bottom [d1,d2], FEP-aware ---
struct DistRestraintTerms {
    int n = 0;
    std::unique_ptr<HostDeviceBuffer<restrdis_t>> recs;  // {ai,aj,ipsi,d1,d2,k}; char labels unused on device
};

// --- restrang: atom triple, harmonic to ang, FEP-aware ---
struct AngRestraintTerms {
    int n = 0;
    std::unique_ptr<HostDeviceBuffer<restrang_t>> recs;  // {ai,aj,ak,ipsi,ang,k}
};
// --- restrwall: group range [ai,aj), Morse/harmonic wall around solvent_center ---
struct WallRestraintTerms {
    int n = 0;
    std::unique_ptr<HostDeviceBuffer<restrwall_t>> recs;  // {ai,aj,d,k,aMorse,dMorse,ih}
};
struct RestraintData {
    PosRestraintTerms posrestr;
    AtomRestraintTerms restrpos;
    SeqRestraintTerms restrseq;
    DistRestraintTerms restrdis;
    AngRestraintTerms restrang;
    WallRestraintTerms restrwall;

    int total() const {
        return posrestr.n + restrpos.n + restrseq.n + restrdis.n + restrang.n + restrwall.n;
    }
    bool enabled() const { return total() > 0; }
};

class RestraintForce {
   public:
    virtual ~RestraintForce() = default;
    void init(Context& ctx, const std::vector<restrpos_t>& restrspos,
              const std::vector<restrseq_t>& restrseqs,
              const std::vector<restrdis_t>& restrdists,
              const std::vector<restrang_t>& restrangs,
              const std::vector<restrwall_t>& restrwalls);
    virtual void calc(Context& ctx) = 0;
    virtual void cleanup() {}
    bool enabled() const { return data_.enabled(); }
    RestraintData& data() { return data_; }

   protected:
    RestraintData data_;
    virtual void init_backend(Context& ctx) {}

   private:
    void build_posrestr(Context& ctx);
    void build_restrpos(Context& ctx, const std::vector<restrpos_t>& restrspos);
    void build_restrseq(Context& ctx, const std::vector<restrseq_t>& restrseqs);
    void build_restrdis(Context& ctx, const std::vector<restrdis_t>& restrdists);
    void build_restrang(Context& ctx, const std::vector<restrang_t>& restrangs);
    void build_restrwall(Context& ctx, const std::vector<restrwall_t>& restrwalls);
};