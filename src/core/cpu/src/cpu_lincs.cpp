#include "cpu_lincs.h"

#include "geometry.h"

void CpuLincs::apply(Context& ctx, HostDeviceBuffer<coord_t>& xcoords) {
    coord_t* candidate = ctx.coords->cpu_data_p;
    const coord_t* reference = xcoords.cpu_data_p;

    apply_to(ctx, candidate, reference);
}

void CpuLincs::init_backend(Context& ctx) {
    int bond_num = data_.n_constraints;
    std::vector<ConstraintBond> constraint_bonds(data_.constraint_bonds->cpu_data_p, data_.constraint_bonds->cpu_data_p + bond_num);
    init_from_bonds(ctx, constraint_bonds);
}

void CpuLincs::init_from_bonds(Context& ctx, const std::vector<ConstraintBond>& bonds) {
    // 1. get all the bond sqrt(inv_wa + inv_wb)
    int bond_num = bonds.size();
    std::vector<double> bond_inv_sqrt_inv_mass_sum(bond_num);
    auto* winv = ctx.winv->cpu_data_p;

    for (int i = 0; i < bond_num; i++) {
        auto& bond = bonds[i];
        int ai = bond.ai - 1, aj = bond.aj - 1;
        bond_inv_sqrt_inv_mass_sum[i] = 1.0 / std::sqrt(winv[ai] + winv[aj]);
    }
    lincs_data_.bond_inv_sqrt_inv_mass_sum = bond_inv_sqrt_inv_mass_sum;

    // 2. get bond graph (which two bond has the same atom)
    std::vector<std::vector<int>> bond_graph(bond_num);
    for (int i = 0; i < bond_num; i++) {
        auto& bond_i = bonds[i];
        int a = bond_i.ai, b = bond_i.aj;
        for (int j = i + 1; j < bond_num; j++) {
            auto& bond_j = bonds[j];
            int a2 = bond_j.ai, b2 = bond_j.aj;
            if (a == a2 || a == b2 || b == a2 || b == b2) {
                bond_graph[i].push_back(j);
                bond_graph[j].push_back(i);
            }
        }
    }
    lincs_data_.bond_graph = bond_graph;

    // 3. calculate the normalized c
    std::vector<std::vector<double>> normalized_c(bond_num);

    auto get_same_atom_in_two_bond = [&](const ConstraintBond& b1, const ConstraintBond& b2) {
        if (b1.ai == b2.ai || b1.ai == b2.aj) {
            return b1.ai;
        } else {
            return b1.aj;
        }
    };

    for (int i = 0; i < bond_num; i++) {
        int sz = bond_graph[i].size();
        if (sz == 0) continue;
        normalized_c[i].resize(sz);
        auto& bond_i = bonds[i];

        for (int j = 0; j < sz; j++) {
            int idx = bond_graph[i][j];
            auto& bond_j = bonds[idx];
            int atom_idx = get_same_atom_in_two_bond(bond_i, bond_j);
            normalized_c[i][j] = -get_sign(atom_idx, bond_i) * get_sign(atom_idx, bond_j) * winv[atom_idx - 1] * bond_inv_sqrt_inv_mass_sum[i] *
                                 bond_inv_sqrt_inv_mass_sum[idx];
        }
    }
    lincs_data_.normalized_c = normalized_c;
}

void CpuLincs::initial_constraint(Context& ctx) {
    /*
    Same idea as cpu_shake.cpp
    */
    auto* coords = ctx.coords->cpu_data_p;
    auto* velocities = ctx.velocities->cpu_data_p;
    std::vector<coord_t> xcoords(ctx.n_atoms);

    for (int i = 0; i < ctx.n_atoms; i++) {
        xcoords[i] = coords[i];
    }

    apply_to(ctx, coords, xcoords.data());

    for (int i = 0; i < ctx.n_atoms; i++) {
        const double dt = ctx.dt;
        xcoords[i].x = coords[i].x - dt * velocities[i].x;
        xcoords[i].y = coords[i].y - dt * velocities[i].y;
        xcoords[i].z = coords[i].z - dt * velocities[i].z;
    }

    apply_to(ctx, xcoords.data(), coords);

    for (int i = 0; i < ctx.n_atoms; i++) {
        const double dt = ctx.dt;
        velocities[i].x = (coords[i].x - xcoords[i].x) / dt;
        velocities[i].y = (coords[i].y - xcoords[i].y) / dt;
        velocities[i].z = (coords[i].z - xcoords[i].z) / dt;
    }
}

int CpuLincs::get_sign(const int a, const ConstraintBond& bond) {
    // we define the direction is bond.ai - bond.aj
    if (a == bond.ai) return 1;
    return -1;
}

void CpuLincs::solve(Context& ctx, std::vector<double> errors, coord_t* coords, const coord_t* xcoords) {
    int bond_num = data_.n_constraints;
    auto& bond_graph = lincs_data_.bond_graph;
    auto& normalized_c = lincs_data_.normalized_c;
    auto& bond_inv_sqrt_inv_mass_sum = lincs_data_.bond_inv_sqrt_inv_mass_sum;
    auto* winv = ctx.winv->cpu_data_p;
    auto* bonds = data_.constraint_bonds->cpu_data_p;

    // now we have C = (I - A)
    std::vector<double> res = errors;
    std::vector<double> current = errors;

    for (int i = 0; i < lincs_settings_.expansion_order; i++) {
        std::vector<double> new_value(bond_num);
        for (int j = 0; j < bond_num; j++) {
            int sz = bond_graph[j].size();
            int a = bonds[j].ai - 1, b = bonds[j].aj - 1;
            coord_t ba = xcoords[a] - xcoords[b];
            ba = ba / norm(ba);

            for (int k = 0; k < sz; k++) {
                int idx = bond_graph[j][k];

                int a2 = bonds[idx].ai - 1, b2 = bonds[idx].aj - 1;
                coord_t bb = xcoords[a2] - xcoords[b2];
                bb = bb / norm(bb);
                new_value[j] += normalized_c[j][k] * current[idx] * dot(ba, bb);
            }
        }
        current = new_value;
        for (int j = 0; j < bond_num; j++) {
            res[j] += current[j];
        }
    }

    // calculate correct lambda
    for (int i = 0; i < bond_num; i++) {
        res[i] = bond_inv_sqrt_inv_mass_sum[i] * res[i];
    }

    // Update coordinates
    for (int i = 0; i < bond_num; i++) {
        auto& bond = bonds[i];
        int ai = bond.ai - 1, aj = bond.aj - 1;
        // old unit bond direction
        coord_t b = (xcoords[ai] - xcoords[aj]);
        b = b / norm(b);

        coord_t dx = -winv[ai] * get_sign(bond.ai, bond) * res[i] * b;

        coords[ai] = coords[ai] + dx;

        dx = -winv[aj] * get_sign(bond.aj, bond) * res[i] * b;
        coords[aj] = coords[aj] + dx;
    }
}

bool CpuLincs::check_ready(const coord_t* coords) {
    int bond_num = data_.n_constraints;
    auto* bonds = data_.constraint_bonds->cpu_data_p;
    for (int i = 0; i < bond_num; i++) {
        auto& bond = bonds[i];
        int ai = bond.ai - 1, aj = bond.aj - 1;

        coord_t cur = coords[ai] - coords[aj];
        double cur_d = norm(cur);
        double aim_d = sqrt(bond.dist2);
        double relative_error = std::fabs(cur_d / aim_d - 1.0);
        if (!std::isfinite(cur_d) || !std::isfinite(relative_error) || relative_error > lincs_settings_.accuracy_tolerance) {
            return false;
        }
    }
    return true;
}

void CpuLincs::apply_to(Context& ctx, coord_t* coords, const coord_t* xcoords) {
    // 1. construct the error array
    int bond_num = data_.n_constraints;
    std::vector<double> errors(bond_num);
    auto* bonds = data_.constraint_bonds->cpu_data_p;
    for (int i = 0; i < bond_num; i++) {
        auto& bond = bonds[i];
        int ai = bond.ai - 1, aj = bond.aj - 1;

        // old unit bond direction
        coord_t b = (xcoords[ai] - xcoords[aj]);
        b = b / norm(b);

        auto ra = coords[ai] - coords[aj];

        double d = sqrt(bonds[i].dist2);  // todo: get the sqrt(dist2) in precompute = =

        errors[i] = dot(b, ra) - d;
    }

    // normalize error
    for (int i = 0; i < bond_num; i++) {
        errors[i] = lincs_data_.bond_inv_sqrt_inv_mass_sum[i] * errors[i];
    }
    solve(ctx, errors, coords, xcoords);

    // Do the correction again
    bool ok = false;
    for (int i = 0; i < lincs_settings_.maximum_rotation_iterations; i++) {
        bool rotation_valid = true;
        for (int j = 0; j < bond_num; j++) {
            int ai = bonds[j].ai - 1, aj = bonds[j].aj - 1;
            coord_t cur = coords[ai] - coords[aj];
            double r2 = norm2(cur);
            double d2 = bonds[j].dist2;
            double radicand = 2.0 * d2 - r2;
            if (!std::isfinite(radicand) || radicand <= 0.0) {
                // LINCS failure
                rotation_valid = false;
                break;
            }
            double p = sqrt(radicand);
            errors[j] = sqrt(d2) - p;

            // normalized error
            errors[j] = lincs_data_.bond_inv_sqrt_inv_mass_sum[j] * errors[j];
        }
        if (!rotation_valid) {
            break;
        }
        solve(ctx, errors, coords, xcoords);
        if (i + 1 >= lincs_settings_.minimum_rotation_iterations && check_ready(coords)) {
            ok = true;
            break;
        }
    }

    if (!ok) {
        for (int i = 0; i < bond_num; i++) {
            int ai = bonds[i].ai - 1, aj = bonds[i].aj - 1;
            coord_t cur = coords[ai] - coords[aj];
            double cur_d = norm(cur);
            double aim_d = sqrt(bonds[i].dist2);
            double relative_error = std::fabs(cur_d / aim_d - 1.0);
            if (!std::isfinite(cur_d) || !std::isfinite(relative_error) || relative_error > lincs_settings_.accuracy_tolerance) {
                printf(">>> Lincs failed, i = %d, j = %d, d = %f, d0 = %f\n", ai, aj, cur_d, aim_d);
            }
        }
        std::exit(EXIT_FAILURE);
    }
}