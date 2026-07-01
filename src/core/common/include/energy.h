#pragma once
#include "host_device_buffer.h"
#include "md_types.h"

enum EnergySlot : int {  // fixed, non-FEP slots
    E_BOND_P_BOND,
    E_BOND_P_ANGLE,
    E_BOND_P_TOR,
    E_BOND_P_IMP,
    E_BOND_W_BOND,
    E_BOND_W_ANGLE,
    E_BOND_W_TOR,
    E_BOND_W_IMP,
    E_NB_PP_COUL,
    E_NB_PP_VDW,
    E_NB_PW_COUL,
    E_NB_PW_VDW,
    E_NB_WW_COUL,
    E_NB_WW_VDW,
    E_RESTR_RADX,
    E_RESTR_POLX,
    E_RESTR_FIX,
    E_RESTR_SHELL,
    E_RESTR_PRES,
    ENERGY_FIXED_COUNT
};

// Per-λ block, appended after the fixed slots (11 values × n_lambdas):
enum EQSlot : int {
    EQ_BOND_BOND,
    EQ_BOND_ANGLE,
    EQ_BOND_TOR,
    EQ_BOND_IMP,
    EQ_NB_QQ_COUL,
    EQ_NB_QQ_VDW,
    EQ_NB_QP_COUL,
    EQ_NB_QP_VDW,
    EQ_NB_QW_COUL,
    EQ_NB_QW_VDW,
    EQ_RESTR_URESTR,
    EQ_STRIDE
};

struct EnergyData {
    // fixed (non-FEP) - fillled directly from slots
    E_bonded_t bond_p, bond_w;
    E_nonbonded_t nb_pp, nb_pw, nb_ww;
    E_restraint_t restraint;  // Uradx, Upolx, Ufix, Ushell, Upres, Urestr

    // per-state (FEP)
    std::vector<E_bonded_t> eq_bond;
    std::vector<E_nonbonded_t> eq_qq, eq_qp, eq_qw;
    std::vector<E_nonbonded_t> eq_qx;  // derived: qq + qp + qw
    std::vector<double> eq_restr;      // slot
    std::vector<double> eq_total;      // derived: per-state Utot

    // derived global / lambda-weighted
    E_bonded_t bond_q;
    E_nonbonded_t nb_qx;
    double Upot = 0, Ukin = 0, Utot = 0;
};

class EnergyBuffer {
    std::unique_ptr<HostDeviceBuffer<energy_accum_t>> accum_;  // length = ENERGY_FIXED_COUNT + EQ_STRIDE * n_lambdas
    EnergyData data_;

    int n_lambdas_;

   public:
    void init(int n_lambdas);
    HD static int eq_index(int n_fixed, int state, EQSlot s)  // slot for (state, field)
    { return n_fixed + state * EQ_STRIDE + s; }
    energy_accum_t* device() { return accum_->gpu_data_p; }  // hand to kernels
    energy_accum_t* host() { return accum_->cpu_data_p; }
    int size() const { return ENERGY_FIXED_COUNT + EQ_STRIDE * n_lambdas_; }

    void reset();
    void download();
    EnergyData& data() { return data_; }
    void unpack();
    void combine(const double* lambdas);
};
