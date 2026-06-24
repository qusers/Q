#include "cpu_nonbonded_force.h"

#include "constants.h"
#include "cpu_force_accumulation.h"

namespace {
void accumulate_energy(Context& ctx, std::vector<energy_accum_t>& e_coul, std::vector<energy_accum_t>& e_vdw, real_t vel, real_t vvdw, uint8_t atom1_type, uint8_t atom2_type, int atom1_state, int atom2_state) {
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

    int n_slots = nb_total_slots(ctx.n_lambdas);
    std::vector<energy_accum_t> e_coul(n_slots), e_vdw(n_slots);

    for (int i = 0; i < sz; i++) {
        const int atom1 = atom_idxs[i];
        if (atom1 == -1) continue;
        const auto& atom1_type = data_.category->cpu_data_p[i];
        const int atom1_state = data_.q_state->cpu_data_p[i];
        const real_t atom1_charge = data_.atom_charge->cpu_data_p[i];
        const vdw_atom_param_t& atom1_vdw = data_.atom_vdw->cpu_data_p[i];
        for (int j = i + 1; j < sz; j++) {
            const int atom2 = atom_idxs[j];
            if (atom2 == -1) continue;
            const auto& atom2_type = data_.category->cpu_data_p[j];
            const int atom2_state = data_.q_state->cpu_data_p[j];
            const auto& bond_type = get_bond_type(ctx.n_atoms_solute, ctx.LJ_matrix->cpu_data_p, atom1, atom1_type, atom2, atom2_type);
            const real_t atom2_charge = data_.atom_charge->cpu_data_p[j];
            const vdw_atom_param_t atom2_vdw = data_.atom_vdw->cpu_data_p[j];

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

            real_t qij = atom1_charge * atom2_charge;
            bool is_14 = (bond_type == BondType::Bond14);
            real_t scaling = is_14 ? ctx.topo.el14_scale : 1;
            real_t2 pair = is_14 ? combine_vdw(ctx.topo.vdw_rule, atom1_vdw.aii_14, atom1_vdw.bii_14, atom2_vdw.aii_14, atom2_vdw.bii_14) : combine_vdw(ctx.topo.vdw_rule, atom1_vdw.aii_normal, atom1_vdw.bii_normal, atom2_vdw.aii_normal, atom2_vdw.bii_normal);

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
