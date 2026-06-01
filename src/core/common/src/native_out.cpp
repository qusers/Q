#include "native_out.h"

#include <algorithm>
#include <cctype>
#include <cstdint>
#include <cstring>
#include <filesystem>
#include <stdexcept>
#include <string>
#include <utility>
#include <vector>

namespace {
std::string trim_copy(const std::string& value) {
    size_t begin = 0;
    while (begin < value.size() && std::isspace(static_cast<unsigned char>(value[begin]))) begin++;

    size_t end = value.size();
    while (end > begin && std::isspace(static_cast<unsigned char>(value[end - 1]))) end--;

    return value.substr(begin, end - begin);
}

std::string lower_copy(std::string value) {
    std::transform(value.begin(), value.end(), value.begin(), [](unsigned char c) {
        return static_cast<char>(std::tolower(c));
    });
    return value;
}

void ensure_parent_dir(const std::string& file_path) {
    std::filesystem::path parent = std::filesystem::path(file_path).parent_path();
    if (!parent.empty()) {
        std::filesystem::create_directories(parent);
    }
}

void write_record(std::ofstream& out, const void* payload, int32_t nbytes) {
    out.write(reinterpret_cast<const char*>(&nbytes), sizeof(nbytes));
    if (nbytes > 0) out.write(reinterpret_cast<const char*>(payload), nbytes);
    out.write(reinterpret_cast<const char*>(&nbytes), sizeof(nbytes));
}

void append_bytes(std::vector<char>& buffer, const void* data, size_t nbytes) {
    const char* bytes = reinterpret_cast<const char*>(data);
    buffer.insert(buffer.end(), bytes, bytes + nbytes);
}

void append_int32(std::vector<char>& buffer, int32_t value) {
    append_bytes(buffer, &value, sizeof(value));
}

void append_double(std::vector<char>& buffer, double value) {
    append_bytes(buffer, &value, sizeof(value));
}

void append_float(std::vector<char>& buffer, float value) {
    append_bytes(buffer, &value, sizeof(value));
}

void append_char_fixed(std::vector<char>& buffer, const std::string& value, size_t width) {
    std::vector<char> field(width, ' ');
    size_t ncopy = std::min(value.size(), width);
    std::memcpy(field.data(), value.data(), ncopy);
    append_bytes(buffer, field.data(), width);
}

void append_bonded(std::vector<char>& record, const E_bonded_t& value) {
    append_double(record, value.Ubond);
    append_double(record, value.Uangle);
    append_double(record, value.Utor);
    append_double(record, value.Uimp);
}

void append_nonbonded(std::vector<char>& record, const E_nonbonded_t& value) {
    append_double(record, value.Ucoul);
    append_double(record, value.Uvdw);
}

void write_restart_record(std::ofstream& out,
                          const coord_t* coords,
                          const vel_t* velocities,
                          int natoms,
                          bool velocity) {
    int32_t nat3 = natoms * 3;
    std::vector<char> payload;
    append_int32(payload, nat3);
    for (int i = 0; i < natoms; i++) {
        double x = velocity ? velocities[i].x : coords[i].x;
        double y = velocity ? velocities[i].y : coords[i].y;
        double z = velocity ? velocities[i].z : coords[i].z;
        append_double(payload, x);
        append_double(payload, y);
        append_double(payload, z);
    }
    write_record(out, payload.data(), static_cast<int32_t>(payload.size()));
}

void write_theta_corr_record(std::ofstream& out, const shell_t* wshells, int n_shells) {
    std::vector<char> payload;
    append_int32(payload, n_shells);
    for (int i = 0; i < n_shells; i++) {
        append_float(payload, static_cast<float>(wshells[i].theta_corr));
    }
    write_record(out, payload.data(), static_cast<int32_t>(payload.size()));
}
}  // namespace

NativeOutput::NativeOutput(NativeOutputConfig config) : config_(std::move(config)) {
    validate_config();
}

void NativeOutput::validate_config() const {
    if (config_.final_file.empty()) {
        throw std::runtime_error("Native output requires a final restart file.");
    }
}

void NativeOutput::validate_runtime_config(Context& ctx) const {
    if (ctx.md.trajectory > 0 && config_.trajectory_file.empty()) {
        throw std::runtime_error("Native output requires a trajectory file when trajectory interval is positive.");
    }
    if (ctx.md.energy > 0 && config_.energy_file.empty()) {
        throw std::runtime_error("Native output requires an energy file when energy interval is positive.");
    }
}

void NativeOutput::ensure_initialized(Context& ctx) {
    if (!output_ready_) {
        init(ctx);
    }
}

void NativeOutput::init(Context& ctx) {
    if (output_ready_) return;

    validate_runtime_config(ctx);
    init_trajectory_mask(ctx);
    open_trajectory(ctx);
    open_energy(ctx);
    output_ready_ = true;
}

void NativeOutput::finish(Context& ctx) {
    ensure_initialized(ctx);
    write_restart_file(ctx);
}

void NativeOutput::shutdown() {
    if (trajectory_stream_.is_open()) trajectory_stream_.close();
    if (energy_stream_.is_open()) energy_stream_.close();
    trajectory_atoms_indices_.clear();
    output_ready_ = false;
}

void NativeOutput::output_trajectory(Context& ctx, int iteration) {
    ensure_initialized(ctx);

    int step = iteration + 1;
    if (ctx.md.trajectory > 0 && step % ctx.md.trajectory == 0) {
        write_trajectory_frame(ctx);
    }
}

void NativeOutput::output_energy(Context& ctx, int iteration) {
    ensure_initialized(ctx);

    int step = iteration + 1;
    if (ctx.md.energy > 0 && step % ctx.md.energy == 0) {
        write_energy_frame(ctx);
    }
}

void NativeOutput::output_restart(Context& ctx, int iteration) {
    ensure_initialized(ctx);

    int step = iteration + 1;
    if (step > 0 && step % 1000 == 0) {
        write_restart_file(ctx);
    }
}

void NativeOutput::init_trajectory_mask(Context& ctx) {
    trajectory_atoms_indices_.clear();
    if (ctx.md.trajectory <= 0) return;

    std::string mask = lower_copy(trim_copy(config_.trajectory_atoms));
    if (mask.empty() || mask == "all") {
        for (int i = 0; i < ctx.n_atoms; i++) {
            trajectory_atoms_indices_.push_back(i);
        }
    } else if (mask == "not excluded") {
        auto* excluded = ctx.excluded->cpu_data_p;
        for (int i = 0; i < ctx.n_atoms; i++) {
            if (!excluded[i]) trajectory_atoms_indices_.push_back(i);
        }
    } else {
        throw std::runtime_error("Unsupported native trajectory_atoms mask: " + config_.trajectory_atoms);
    }
}

void NativeOutput::open_trajectory(Context& ctx) {
    if (ctx.md.trajectory <= 0) return;

    ensure_parent_dir(config_.trajectory_file);
    trajectory_stream_.open(config_.trajectory_file.c_str(), std::ios::binary | std::ios::trunc);
    if (!trajectory_stream_) {
        throw std::runtime_error("Could not write native trajectory file " + config_.trajectory_file);
    }

    std::vector<char> rec1;
    append_char_fixed(rec1, "CORD", 4);
    append_int32(rec1, ctx.md.steps / ctx.md.trajectory);
    append_int32(rec1, ctx.md.trajectory);
    append_int32(rec1, ctx.md.trajectory);
    append_int32(rec1, ctx.md.steps);
    append_int32(rec1, 0);
    append_int32(rec1, 0);
    append_int32(rec1, 0);
    append_int32(rec1, static_cast<int32_t>(ctx.Ndegfree));
    append_int32(rec1, 0);
    append_int32(rec1, 0);
    append_int32(rec1, 0);
    append_int32(rec1, 0);
    append_int32(rec1, 0);
    append_int32(rec1, 0);
    append_int32(rec1, 0);
    append_int32(rec1, 0);
    append_int32(rec1, 0);
    append_int32(rec1, 0);
    append_int32(rec1, 0);
    append_int32(rec1, 0);
    write_record(trajectory_stream_, rec1.data(), static_cast<int32_t>(rec1.size()));

    std::string mask_line = config_.trajectory_atoms.empty() ? "all" : config_.trajectory_atoms;
    int32_t mask_rows = 1;
    std::vector<char> rec2;
    append_int32(rec2, 2 + mask_rows);
    append_char_fixed(rec2, "Q DCD trajectory version 6.0.1", 80);
    append_char_fixed(rec2, config_.topology_file, 80);
    append_char_fixed(rec2, mask_line, 80);
    write_record(trajectory_stream_, rec2.data(), static_cast<int32_t>(rec2.size()));

    int32_t included = static_cast<int32_t>(trajectory_atoms_indices_.size());
    write_record(trajectory_stream_, &included, sizeof(included));
}

void NativeOutput::open_energy(Context& ctx) {
    if (ctx.md.energy <= 0) return;

    ensure_parent_dir(config_.energy_file);
    energy_stream_.open(config_.energy_file.c_str(), std::ios::binary | std::ios::trunc);
    if (!energy_stream_) {
        throw std::runtime_error("Could not write native energy file " + config_.energy_file);
    }
}

void NativeOutput::write_trajectory_frame(Context& ctx) {
    if (ctx.md.trajectory <= 0 || !trajectory_stream_) return;

    std::vector<float> axis(trajectory_atoms_indices_.size());
    auto* coords = ctx.coords->cpu_data_p;
    for (size_t i = 0; i < trajectory_atoms_indices_.size(); i++) axis[i] = static_cast<float>(coords[trajectory_atoms_indices_[i]].x);
    write_record(trajectory_stream_, axis.data(), static_cast<int32_t>(axis.size() * sizeof(float)));
    for (size_t i = 0; i < trajectory_atoms_indices_.size(); i++) axis[i] = static_cast<float>(coords[trajectory_atoms_indices_[i]].y);
    write_record(trajectory_stream_, axis.data(), static_cast<int32_t>(axis.size() * sizeof(float)));
    for (size_t i = 0; i < trajectory_atoms_indices_.size(); i++) axis[i] = static_cast<float>(coords[trajectory_atoms_indices_[i]].z);
    write_record(trajectory_stream_, axis.data(), static_cast<int32_t>(axis.size() * sizeof(float)));
}

void NativeOutput::write_energy_frame(Context& ctx) {
    if (ctx.md.energy <= 0 || !energy_stream_) return;

    auto* lambdas = ctx.lambdas ? ctx.lambdas->cpu_data_p : nullptr;
    auto* EQ_restraint = ctx.EQ_restraint->cpu_data_p;
    for (int state = 0; state < ctx.n_lambdas; state++) {
        std::vector<char> record;
        append_int32(record, state + 1);
        append_double(record, lambdas[state]);
        append_double(record, ctx.EQ_total[state].Utot);
        append_bonded(record, ctx.EQ_bond[state]);
        append_nonbonded(record, ctx.EQ_nonbond_qx[state]);
        append_nonbonded(record, ctx.EQ_nonbond_qq[state]);
        append_nonbonded(record, ctx.EQ_nonbond_qp[state]);
        append_nonbonded(record, ctx.EQ_nonbond_qw[state]);
        append_double(record, EQ_restraint[state].Urestr);
        write_record(energy_stream_, record.data(), static_cast<int32_t>(record.size()));
    }

    write_record(energy_stream_, nullptr, 0);
}

void NativeOutput::write_restart_file(Context& ctx) const {
    ensure_parent_dir(config_.final_file);
    std::ofstream out(config_.final_file.c_str(), std::ios::binary | std::ios::trunc);
    if (!out) {
        throw std::runtime_error("Could not write native final restart file " + config_.final_file);
    }

    write_restart_record(out, ctx.coords->cpu_data_p, nullptr, ctx.n_atoms, false);
    write_restart_record(out, nullptr, ctx.velocities->cpu_data_p, ctx.n_atoms, true);
    if (ctx.md.polarisation && ctx.wshells && ctx.n_shells > 0) {
        write_theta_corr_record(out, ctx.wshells->cpu_data_p, ctx.n_shells);
    }
}
