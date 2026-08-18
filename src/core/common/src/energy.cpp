#include "energy.h"

#include "cpu_force_accumulation.h"

void EnergyBuffer::init(int n_lambdas) {
    n_lambdas_ = n_lambdas;
    accum_ = std::make_unique<HostDeviceBuffer<energy_accum_t>>(size());
    data_.eq_bond.resize(n_lambdas_);
    data_.eq_qq.resize(n_lambdas_);
    data_.eq_qp.resize(n_lambdas_);
    data_.eq_qw.resize(n_lambdas_);
    data_.eq_qx.resize(n_lambdas_);
    data_.eq_restr.resize(n_lambdas_);
    data_.eq_total.resize(n_lambdas_);
}

void EnergyBuffer::reset() {
    accum_->zero();
}

void EnergyBuffer::download() {
    accum_->download();
}
// Decode the flat fixed-point slots into named double fields. This only does
// work that does NOT depend on lambda: raw per-category energies and the
// per-state derived values (eq_qx = qq+qp+qw, eq_total = per-state sum).
// The lambda-weighted globals (bond_q, nb_qx, restraint.Urestr, Upot/Utot)
// and the lambda==0 zeroing are finalized separately, where lambda is known.
void EnergyBuffer::unpack() {
    const energy_accum_t* a = host();
    auto E = [&](int slot) { return energy_from_accum(a[slot]); };
    auto EQ = [&](int state, EQSlot s) { return E(eq_index(ENERGY_FIXED_COUNT, state, s)); };
    auto& d = data_;

    // fixed bonded
    d.bond_p = {E(E_BOND_P_BOND), E(E_BOND_P_ANGLE), E(E_BOND_P_TOR), E(E_BOND_P_IMP)};
    d.bond_w = {E(E_BOND_W_BOND), E(E_BOND_W_ANGLE), E(E_BOND_W_TOR), E(E_BOND_W_IMP)};

    // fixed nonbonded
    d.nb_pp = {E(E_NB_PP_COUL), E(E_NB_PP_VDW)};
    d.nb_pw = {E(E_NB_PW_COUL), E(E_NB_PW_VDW)};
    d.nb_ww = {E(E_NB_WW_COUL), E(E_NB_WW_VDW)};

    // fixed restraint components (raw). Upres here is only the direct kernel
    // contribution; the lambda-weighted per-state term is added in combine.
    d.restraint.Uradx = E(E_RESTR_RADX);
    d.restraint.Upolx = E(E_RESTR_POLX);
    d.restraint.Ufix = E(E_RESTR_FIX);
    d.restraint.Ushell = E(E_RESTR_SHELL);
    d.restraint.Upres = E(E_RESTR_PRES);

    // per-state (unweighted)
    for (int s = 0; s < n_lambdas_; s++) {
        d.eq_bond[s] = {EQ(s, EQ_BOND_BOND), EQ(s, EQ_BOND_ANGLE), EQ(s, EQ_BOND_TOR), EQ(s, EQ_BOND_IMP)};
        d.eq_qq[s] = {EQ(s, EQ_NB_QQ_COUL), EQ(s, EQ_NB_QQ_VDW)};
        d.eq_qp[s] = {EQ(s, EQ_NB_QP_COUL), EQ(s, EQ_NB_QP_VDW)};
        d.eq_qw[s] = {EQ(s, EQ_NB_QW_COUL), EQ(s, EQ_NB_QW_VDW)};
        d.eq_restr[s] = EQ(s, EQ_RESTR_URESTR);

        d.eq_qx[s] = {d.eq_qq[s].Ucoul + d.eq_qp[s].Ucoul + d.eq_qw[s].Ucoul,
                      d.eq_qq[s].Uvdw + d.eq_qp[s].Uvdw + d.eq_qw[s].Uvdw};
        d.eq_total[s] = d.eq_bond[s].Ubond + d.eq_bond[s].Uangle + d.eq_bond[s].Utor + d.eq_bond[s].Uimp +
                        d.eq_qx[s].Ucoul + d.eq_qx[s].Uvdw + d.eq_restr[s];
    }
}

void EnergyBuffer::combine(const double* lambdas) {
    auto& d = data_;
    d.bond_q = {};
    d.nb_qx = {};
    for (int s = 0; s < n_lambdas_; s++) {
        double L = lambdas[s];
        if (L == 0.0) {  // zero-weight state contributes nothing
            continue;
        }
        d.bond_q.Ubond += d.eq_bond[s].Ubond * L;
        d.bond_q.Uangle += d.eq_bond[s].Uangle * L;
        d.bond_q.Utor += d.eq_bond[s].Utor * L;
        d.bond_q.Uimp += d.eq_bond[s].Uimp * L;
        d.nb_qx.Ucoul += d.eq_qx[s].Ucoul * L;
        d.nb_qx.Uvdw += d.eq_qx[s].Uvdw * L;
        d.restraint.Upres += d.eq_restr[s] * L;  // λ-weighted per-state restraint
    }
    d.restraint.Urestr = d.restraint.Uradx + d.restraint.Upolx +
                         d.restraint.Ushell + d.restraint.Ufix + d.restraint.Upres;
    d.Upot = d.bond_p.Ubond + d.bond_w.Ubond + d.bond_p.Uangle + d.bond_w.Uangle +
             d.bond_p.Utor + d.bond_p.Uimp +
             d.nb_pp.Ucoul + d.nb_pp.Uvdw + d.nb_pw.Ucoul + d.nb_pw.Uvdw +
             d.nb_ww.Ucoul + d.nb_ww.Uvdw +
             d.bond_q.Ubond + d.bond_q.Uangle + d.bond_q.Utor + d.bond_q.Uimp +
             d.nb_qx.Ucoul + d.nb_qx.Uvdw + d.restraint.Urestr;
    d.Utot = d.Upot + d.Ukin;
}