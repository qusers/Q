#include "q_output.h"
#include "q_input.h"

#include "context.h"

#include <algorithm>
#include <cctype>
#include <cerrno>
#include <cstdint>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <fstream>
#include <string>
#include <sys/stat.h>
#include <vector>

static std::string q_final_file;
static std::string q_trajectory_file;
static std::string q_energy_file;
static std::string q_trajectory_atoms;
static std::string q_topology_file;
static int q_steps = 0;
static int q_trajectory_interval = 0;
static int q_energy_interval = 0;
static std::ofstream q_traj_stream;
static std::ofstream q_energy_stream;
static std::vector<int> q_traj_atoms;
static bool q_output_ready = false;

static void fatal(const std::string &message) {
    printf(">>> FATAL: %s\n", message.c_str());
    exit(EXIT_FAILURE);
}

static std::string dirname_of(const std::string &path) {
    size_t slash = path.find_last_of('/');
    if (slash == std::string::npos) return ".";
    if (slash == 0) return "/";
    return path.substr(0, slash);
}

static void make_dir_recursive(const std::string &path) {
    if (path.empty() || path == ".") return;
    std::string current;
    size_t pos = 0;
    if (path[0] == '/') {
        current = "/";
        pos = 1;
    }
    while (pos <= path.size()) {
        size_t next = path.find('/', pos);
        std::string part = path.substr(pos, next == std::string::npos ? std::string::npos : next - pos);
        if (!part.empty()) {
            if (current.size() > 1 && current[current.size() - 1] != '/') current += "/";
            current += part;
            if (mkdir(current.c_str(), 0775) != 0 && errno != EEXIST) {
                fatal("Could not create directory " + current);
            }
        }
        if (next == std::string::npos) break;
        pos = next + 1;
    }
}

static std::string trim(std::string value) {
    size_t first = value.find_first_not_of(" \t\r\n");
    if (first == std::string::npos) return "";
    size_t last = value.find_last_not_of(" \t\r\n");
    return value.substr(first, last - first + 1);
}

static std::string lower(std::string value) {
    std::transform(value.begin(), value.end(), value.begin(), [](unsigned char c) {
        return (char) std::tolower(c);
    });
    return value;
}

static void write_record(std::ofstream &out, const void *payload, int32_t nbytes) {
    out.write((const char*) &nbytes, sizeof(nbytes));
    if (nbytes > 0) out.write((const char*) payload, nbytes);
    out.write((const char*) &nbytes, sizeof(nbytes));
}

static void append_bytes(std::vector<char> &buffer, const void *data, size_t nbytes) {
    const char *bytes = (const char*) data;
    buffer.insert(buffer.end(), bytes, bytes + nbytes);
}

static void append_int32(std::vector<char> &buffer, int32_t value) {
    append_bytes(buffer, &value, sizeof(value));
}

static void append_double(std::vector<char> &buffer, double value) {
    append_bytes(buffer, &value, sizeof(value));
}

static void append_char_fixed(std::vector<char> &buffer, const std::string &value, size_t width) {
    std::vector<char> field(width, ' ');
    size_t ncopy = value.size() < width ? value.size() : width;
    memcpy(field.data(), value.data(), ncopy);
    append_bytes(buffer, field.data(), width);
}

static void write_restart_record(std::ofstream &out, const coord_t *coords_in, const vel_t *vels_in, int natoms, bool velocity) {
    int32_t nat3 = natoms * 3;
    std::vector<char> payload;
    append_int32(payload, nat3);
    for (int i = 0; i < natoms; i++) {
        double x = velocity ? vels_in[i].x : coords_in[i].x;
        double y = velocity ? vels_in[i].y : coords_in[i].y;
        double z = velocity ? vels_in[i].z : coords_in[i].z;
        append_double(payload, x);
        append_double(payload, y);
        append_double(payload, z);
    }
    write_record(out, payload.data(), (int32_t) payload.size());
}

void q_output_write_restart(const coord_t *final_coords, const vel_t *final_velocities, int natoms) {
    if (q_input_mode != Q_INPUT_NATIVE || q_final_file.empty()) return;
    make_dir_recursive(dirname_of(q_final_file));
    std::ofstream out(q_final_file.c_str(), std::ios::binary | std::ios::trunc);
    if (!out) fatal("Could not write final restart file " + q_final_file);
    write_restart_record(out, final_coords, NULL, natoms, false);
    write_restart_record(out, NULL, final_velocities, natoms, true);
}

static void init_trajectory_mask() {
    auto& ctx = Context::instance();
    q_traj_atoms.clear();
    if (q_trajectory_interval <= 0) return;

    std::string mask = lower(trim(q_trajectory_atoms));
    if (mask.empty() || mask == "all") {
        for (int i = 0; i < ctx.n_atoms; i++) q_traj_atoms.push_back(i);
    }
    else if (mask == "not excluded") {
        auto *excluded = ctx.excluded->cpu_data_p;
        for (int i = 0; i < ctx.n_atoms; i++) {
            if (!excluded[i]) q_traj_atoms.push_back(i);
        }
    }
    else {
        fatal("Unsupported trajectory_atoms mask for native QGPU output: " + q_trajectory_atoms);
    }
}

static void open_trajectory() {
    auto& ctx = Context::instance();
    if (q_trajectory_interval <= 0) return;
    if (q_trajectory_file.empty()) fatal("[files] trajectory is required when trajectory interval is positive.");

    make_dir_recursive(dirname_of(q_trajectory_file));
    q_traj_stream.open(q_trajectory_file.c_str(), std::ios::binary | std::ios::trunc);
    if (!q_traj_stream) fatal("Could not write trajectory file " + q_trajectory_file);

    std::vector<char> rec1;
    append_char_fixed(rec1, "CORD", 4);
    append_int32(rec1, q_steps / q_trajectory_interval);
    append_int32(rec1, q_trajectory_interval);
    append_int32(rec1, q_trajectory_interval);
    append_int32(rec1, q_steps);
    append_int32(rec1, 0);
    append_int32(rec1, 0);
    append_int32(rec1, 0);
    append_int32(rec1, (int32_t) ctx.Ndegfree);
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
    write_record(q_traj_stream, rec1.data(), (int32_t) rec1.size());

    std::string mask_line = q_trajectory_atoms.empty() ? "all" : q_trajectory_atoms;
    int32_t mask_rows = 1;
    std::vector<char> rec2;
    append_int32(rec2, 2 + mask_rows);
    append_char_fixed(rec2, "Q DCD trajectory version 6.0.1", 80);
    append_char_fixed(rec2, q_topology_file, 80);
    append_char_fixed(rec2, mask_line, 80);
    write_record(q_traj_stream, rec2.data(), (int32_t) rec2.size());

    int32_t included = (int32_t) q_traj_atoms.size();
    write_record(q_traj_stream, &included, sizeof(included));
}

static void write_trajectory_frame() {
    auto& ctx = Context::instance();
    if (q_trajectory_interval <= 0 || !q_traj_stream) return;

    std::vector<float> axis(q_traj_atoms.size());
    auto *coords = ctx.coords->cpu_data_p;
    for (size_t i = 0; i < q_traj_atoms.size(); i++) axis[i] = (float) coords[q_traj_atoms[i]].x;
    write_record(q_traj_stream, axis.data(), (int32_t) (axis.size() * sizeof(float)));
    for (size_t i = 0; i < q_traj_atoms.size(); i++) axis[i] = (float) coords[q_traj_atoms[i]].y;
    write_record(q_traj_stream, axis.data(), (int32_t) (axis.size() * sizeof(float)));
    for (size_t i = 0; i < q_traj_atoms.size(); i++) axis[i] = (float) coords[q_traj_atoms[i]].z;
    write_record(q_traj_stream, axis.data(), (int32_t) (axis.size() * sizeof(float)));
}

static void open_energy() {
    if (q_energy_interval <= 0) return;
    if (q_energy_file.empty()) fatal("[files] energy is required when energy interval is positive.");

    make_dir_recursive(dirname_of(q_energy_file));
    q_energy_stream.open(q_energy_file.c_str(), std::ios::binary | std::ios::trunc);
    if (!q_energy_stream) fatal("Could not write energy file " + q_energy_file);
}

static void append_bonded(std::vector<char> &record, const E_bonded_t &value) {
    append_double(record, value.Ubond);
    append_double(record, value.Uangle);
    append_double(record, value.Utor);
    append_double(record, value.Uimp);
}

static void append_nonbonded(std::vector<char> &record, const E_nonbonded_t &value) {
    append_double(record, value.Ucoul);
    append_double(record, value.Uvdw);
}

static void write_energy_frame() {
    auto& ctx = Context::instance();
    if (q_energy_interval <= 0 || !q_energy_stream) return;

    auto *lambdas = ctx.lambdas ? ctx.lambdas->cpu_data_p : nullptr;
    auto *EQ_restraint = ctx.EQ_restraint->cpu_data_p;
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
        write_record(q_energy_stream, record.data(), (int32_t) record.size());
    }

    write_record(q_energy_stream, NULL, 0);
}

void q_output_configure(const char *final_file,
                        const char *trajectory_file,
                        const char *energy_file,
                        const char *trajectory_atoms,
                        const char *topology_file,
                        int steps,
                        int trajectory_interval,
                        int energy_interval) {
    q_final_file = final_file ? final_file : "";
    q_trajectory_file = trajectory_file ? trajectory_file : "";
    q_energy_file = energy_file ? energy_file : "";
    q_trajectory_atoms = trajectory_atoms ? trajectory_atoms : "";
    q_topology_file = topology_file ? topology_file : "";
    q_steps = steps;
    q_trajectory_interval = trajectory_interval;
    q_energy_interval = energy_interval;

    if (q_final_file.empty()) fatal("Q input [files] section must specify final.");
    if (q_trajectory_interval > 0 && q_trajectory_file.empty()) fatal("Q input [files] section must specify trajectory when trajectory interval is positive.");
    if (q_energy_interval > 0 && q_energy_file.empty()) fatal("Q input [files] section must specify energy when energy interval is positive.");
}

void q_output_init() {
    if (q_input_mode != Q_INPUT_NATIVE || q_output_ready) return;
    init_trajectory_mask();
    open_trajectory();
    open_energy();
    q_output_ready = true;
}

void q_output_step(int iteration) {
    auto& ctx = Context::instance();
    if (q_input_mode != Q_INPUT_NATIVE) return;
    if (!q_output_ready) q_output_init();

    int step = iteration + 1;
    if (q_trajectory_interval > 0 && step % q_trajectory_interval == 0) {
        write_trajectory_frame();
    }
    if (q_energy_interval > 0 && step % q_energy_interval == 0) {
        write_energy_frame();
    }
    if (step > 0 && step % 1000 == 0) {
        q_output_write_restart(ctx.coords->cpu_data_p, ctx.velocities->cpu_data_p, ctx.n_atoms);
    }
}

void q_output_finish() {
    auto& ctx = Context::instance();
    if (q_input_mode != Q_INPUT_NATIVE) return;
    if (!q_output_ready) q_output_init();
    q_output_write_restart(ctx.coords->cpu_data_p, ctx.velocities->cpu_data_p, ctx.n_atoms);
}

void q_output_shutdown() {
    if (q_traj_stream.is_open()) q_traj_stream.close();
    if (q_energy_stream.is_open()) q_energy_stream.close();
    q_traj_atoms.clear();
    q_output_ready = false;
}
