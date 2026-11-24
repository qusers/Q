#pragma once
#include "host_device_buffer.h"
#include "md_types.h"
// Foreign-work targets per stage ([foreign_lambdas]); fixed so that the energy buffer keeps its size.
constexpr int MAX_FOREIGN_LAMBDAS = 8;

// Passed by value to the pair kernels. n == 0 means "do not evaluate", which is every step
// that does not write energies.
struct foreign_targets_t {
    int n = 0;
    double weight[MAX_FOREIGN_LAMBDAS][2] = {};  // lambda' of state 1 and 2
    double sc[MAX_FOREIGN_LAMBDAS][2] = {};      // (1 - lambda')^(1/6)

    double rest = 1.0;
    double rest_sqrt = 1.0;
    double rest_target[MAX_FOREIGN_LAMBDAS] = {};
    double rest_target_sqrt[MAX_FOREIGN_LAMBDAS] = {};
};

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

    E_LRF,

    E_RESTR_RADX,
    E_RESTR_POLX,
    E_RESTR_FIX,
    E_RESTR_SHELL,
    E_RESTR_PRES,

    E_BOND_REST_TOR,  // topology torsions with four Q atoms; raw, scaled by rest in combine()

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
    double lrf = 0;
    double rest_tor = 0;  // raw energy of the torsions scaled by REST (E_BOND_REST_TOR)

    // per-state (FEP)
    std::vector<E_bonded_t> eq_bond;
    std::vector<E_nonbonded_t> eq_qq, eq_qp, eq_qw;
    std::vector<E_nonbonded_t> eq_qx;  // derived: qq + qp + qw
    std::vector<double> eq_restr;      // slot
    std::vector<double> eq_total;      // derived: per-state Utot

    // derived global / lambda-weighted
    E_bonded_t bond_q;
    E_nonbonded_t nb_qx;

    // Per target: sum over states of lambda'_s * sum over soft pairs of [V_s(lambda'_s) - V_s(lambda_s)].
    std::vector<double> foreign_shape;
    // Per target: U(lambda') - U(lambda) at the current coordinates.
    std::vector<double> foreign_work;

    double Upot = 0, Ukin = 0, Utot = 0;
};

class EnergyBuffer {
    std::unique_ptr<HostDeviceBuffer<energy_accum_t>> accum_;  // length = ENERGY_FIXED_COUNT + EQ_STRIDE * n_lambdas + MAX_FOREIGN_LAMBDAS
    EnergyData data_;

    int n_lambdas_;

   public:
    void init(int n_lambdas);
    HD static int eq_index(int n_fixed, int state, EQSlot s)  // slot for (state, field)
    { return n_fixed + state * EQ_STRIDE + s; }
    energy_accum_t* device() { return accum_->gpu_data_p; }  // hand to kernels
    energy_accum_t* host() { return accum_->cpu_data_p; }
    int size() const { return ENERGY_FIXED_COUNT + EQ_STRIDE * n_lambdas_ + MAX_FOREIGN_LAMBDAS; }
    HD static int foreign_index(int n_lambdas, int k) { return ENERGY_FIXED_COUNT + n_lambdas * EQ_STRIDE + k; }

    void reset();
    void download();
    EnergyData& data() { return data_; }
    void unpack();
    void combine(const double* lambdas, double rest = 1.0);
};
