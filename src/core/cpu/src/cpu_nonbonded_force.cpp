#include "cpu_nonbonded_force.h"

#include "constants.h"
#include "cpu_force_accumulation.h"

namespace {
enum class BondType : uint8_t { Bond23,
                                Bond14,
                                NonBond };

BondType get_bond_type(Context& ctx, int atom1, uint8_t atom1_type, int atom2, uint8_t atom2_type) {
    bool w1 = atom1_type == static_cast<uint8_t>(AtomCategory::W);
    bool w2 = atom2_type == static_cast<uint8_t>(AtomCategory::W);

    if (w1 != w2) {
        return BondType::NonBond;  // solute-water never bonded
    }

    if (w1) {
        // both w
        if ((atom1 - ctx.n_atoms_solute) / 3 == (atom2 - ctx.n_atoms_solute) / 3) {
            return BondType::Bond23;
        }
        return BondType::NonBond;
    }

    int lj_matrix_value = ctx.LJ_matrix->cpu_data_p[atom1 * ctx.n_atoms_solute + atom2];
    if (lj_matrix_value == 3) {
        return BondType::Bond23;
    } else if (lj_matrix_value == 1) {
        return BondType::Bond14;
    }
    return BondType::NonBond;
}

std::pair<real_t, real_t> calc_electrostatic(real_t qij, real_t coulomb_constant, real_t inv_dis) {
    if (qij == 0) return {0, 0};
    real_t vel = qij * coulomb_constant * inv_dis;  // k * qi * qj / r
    real_t dvel = -vel * inv_dis;                   /// -k * qi * qj / r2
    return {vel, dvel};
}

std::pair<real_t, real_t> calc_vdw(const std::pair<real_t, real_t> pair, real_t inv_dis) {
    if (pair.first == 0 && pair.second == 0) return {0, 0};
    real_t inv_dis3 = inv_dis * inv_dis * inv_dis;
    real_t inv_dis6 = inv_dis3 * inv_dis3;

    real_t V_a = pair.first * inv_dis6 * inv_dis6;  // precomputed A / r^12
    real_t V_b = pair.second * inv_dis6;            // precomputed B / r^6

    real_t vvdw = V_a - V_b;
    real_t dvvdw = (static_cast<real_t>(-12.0) * V_a + static_cast<real_t>(6.0) * V_b) * inv_dis;

    return {vvdw, dvvdw};
}

void accumulate_energy(Context &ctx, std::vector<energy_accum_t> &e_coul, std::vector<energy_accum_t> &e_vdw, real_t vel, real_t vvdw, uint8_t atom1_type, uint8_t atom2_type, int atom1_state, int atom2_state) {
    int slot = nb_energy_slot(atom1_type, atom2_type, atom1_state, atom2_state, ctx.n_lambdas); 
    add_energy(e_coul[slot], vel);
    add_energy(e_vdw[slot], vvdw);
}

}  // namespace

void CpuNonbondedForce::calc(Context& ctx) {
    const auto& atom_idxs = data_.atom_idx->cpu_data_p;
    const auto& coords = ctx.coords->cpu_data_p;
    auto& dvelocities = ctx.dvelocities->cpu_data_p;
    int sz = data_.n_total;
    int n_charge_type = data_.n_charge_types;
    int n_catype_type = data_.n_catype_types;

    int n_slots = nb_total_slots(ctx.n_lambdas);
    std::vector<energy_accum_t> e_coul(n_slots), e_vdw(n_slots);

    for (int i = 0; i < sz; i++) {
        const int atom1 = atom_idxs[i];
        const auto& atom1_type = data_.category->cpu_data_p[i];
        const int atom1_state = data_.q_state->cpu_data_p[i];
        const int atom1_charge_type = data_.charge_types->cpu_data_p[i];
        const int atom1_catype = data_.catype_types->cpu_data_p[i];
        for (int j = i + 1; j < sz; j++) {
            const int atom2 = atom_idxs[j];
            const auto& atom2_type = data_.category->cpu_data_p[j];
            const int atom2_state = data_.q_state->cpu_data_p[j];
            const auto& bond_type = get_bond_type(ctx, atom1, atom1_type, atom2, atom2_type);
            const int atom2_charge_type = data_.charge_types->cpu_data_p[j];
            const int atom2_catype = data_.catype_types->cpu_data_p[j];

            if (bond_type == BondType::Bond23) continue;
            if (atom1_type == static_cast<uint8_t>(AtomCategory::Q) && atom2_type == static_cast<uint8_t>(AtomCategory::Q) && atom1_state != atom2_state) {
                continue;
            }

            real_t dx = coords[atom2].x - coords[atom1].x;
            real_t dy = coords[atom2].y - coords[atom1].y;
            real_t dz = coords[atom2].z - coords[atom1].z;
            real_t dis2 = dx * dx + dy * dy + dz * dz;
            real_t inv_dis2 = static_cast<real_t>(1.0) / dis2;
            real_t inv_dis = sqrt(inv_dis2);

            real_t qij = data_.charge_pair_products->cpu_data_p[n_charge_type * atom1_charge_type + atom2_charge_type];
            vdw_pair_param_t vdw_pair = data_.catype_pair_params->cpu_data_p[n_catype_type * atom1_catype + atom2_catype];
            real_t scaling = bond_type == BondType::Bond14 ? ctx.topo.el14_scale : 1;
            std::pair<real_t, real_t> pair = (bond_type == BondType::Bond14) ? std::make_pair(vdw_pair.a_14, vdw_pair.b_14) : std::make_pair(vdw_pair.a_normal, vdw_pair.b_normal);

            auto [vel, dvel] = calc_electrostatic(qij * scaling, ctx.topo.coulomb_constant, inv_dis);
            auto [vvdw, dvvdw] = calc_vdw(pair, inv_dis);

            real_t lambda = std::min(data_.atom_lambdas->cpu_data_p[i], data_.atom_lambdas->cpu_data_p[j]);

            real_t dva = (dvel + dvvdw) * inv_dis * lambda;

            add_force(dvelocities[atom1].x, -dva * dx);
            add_force(dvelocities[atom1].y, -dva * dy);
            add_force(dvelocities[atom1].z, -dva * dz);

            add_force(dvelocities[atom2].x, dva * dx);
            add_force(dvelocities[atom2].y, dva * dy);
            add_force(dvelocities[atom2].z, dva * dz);

            // Accumulate energy
            accumulate_energy(ctx, e_coul, e_vdw, vel, vvdw, atom1_type, atom2_type, atom1_state, atom2_state);
        }
    }

    auto store = [&](E_nonbonded_t& e, int slot) {
        e.Ucoul = energy_from_accum(e_coul[slot]);
        e.Uvdw = energy_from_accum(e_vdw[slot]);
    };

    store(ctx.E_nonbond_pp, NB_PP);
    store(ctx.E_nonbond_pw, NB_PW);
    store(ctx.E_nonbond_ww, NB_WW);

    for (int s = 0; s < ctx.n_lambdas; s++) {
        store(ctx.EQ_nonbond_qq[s], nb_qq_slot(s, ctx.n_lambdas));
        store(ctx.EQ_nonbond_qp[s], nb_qp_slot(s, ctx.n_lambdas));
        store(ctx.EQ_nonbond_qw[s], nb_qw_slot(s, ctx.n_lambdas));
    }
}
