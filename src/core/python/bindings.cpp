#include <cuda_runtime.h>
#include <pybind11/numpy.h>
#include <pybind11/pybind11.h>
#include <pybind11/stl.h>

#include <filesystem>
#include <memory>
#include <string>
#include <vector>

#include "command_parser.h"
#include "cuda_handler.cuh"
#include "cuda_runtime_utility.h"
#include "engine_types.h"
#include "inp_parser.h"
#include "qgpu_error.h"

namespace py = pybind11;

namespace {

void select_device(int device) {
    int count = 0;
    check_cuda(cudaGetDeviceCount(&count));
    if (device < 0 || device >= count)
        throw InputError("Invalid visible CUDA device index: " + std::to_string(device));
    check_cuda(cudaSetDevice(device));
}

std::string abs_inp(const std::string& file) {
    const auto p = std::filesystem::absolute(file);
    if (!std::filesystem::is_regular_file(p)) throw InputError("Input file does not exist: " + p.string());
    if (p.extension() != ".inp") throw InputError("Expected a Q .inp input file.");
    return p.string();
}

// Parsed settings as qgpu.Config field values (the parser's defaults included).
// Restraints are value tuples in the argument order of the qgpu.restraints classes.
py::dict settings_dict(const ParseResult& p) {
    const md_t& md = p.md;
    const NativeOutputConfig& out = p.native_output;

    py::dict d;
    d["seed"] = md.random_seed;
    // None: not given, so no new velocities (a restart file is required)
    d["initial_temperature"] = md.generate_velocities ? py::object(py::float_(md.initial_temperature)) : py::object(py::none());
    d["shell_radius"] = md.shell_radius;
    d["exclude_bonded"] = md.exclude_bonded;
    d["solute_solute"] = md.solute_solute;
    d["solvent_solvent"] = md.solvent_solvent;
    d["solute_solvent"] = md.solute_solvent;
    d["q_atom"] = md.q_atom;
    d["lrf"] = md.lrf;
    d["lrf_cutoff"] = md.lrf_cutoff;
    d["radial_force"] = md.radial_force;
    d["polarisation"] = md.polarisation;
    d["polarisation_force"] = md.polarisation_force;
    d["charge_correction"] = md.charge_correction;
    d["non_bond"] = md.non_bond;
    d["hmr"] = md.hmr;
    d["hmr_target_mass"] = md.hmr_target_mass;
    d["solute_constraint_algorithm"] = md.solute_constraint_algorithm;
    d["solvent_constraint_algorithm"] = md.solvent_constraint_algorithm;

    d["steps"] = md.steps;
    d["stepsize"] = md.stepsize;
    d["temperature"] = md.temperature;
    d["thermostat"] = md.thermostat;
    d["bath_coupling"] = md.bath_coupling;
    d["constrain_solvent"] = md.constraint_solvent;
    d["constrain_solute"] = md.constraint_solute;
    d["constrain_hydrogens"] = md.constraint_hydrogens;
    d["separate_scaling"] = md.separate_scaling;
    d["shell_force"] = md.shell_force;
    d["output"] = md.output;
    d["energy"] = md.energy;
    d["trajectory"] = md.trajectory;
    d["energy_file"] = out.energy_file;
    d["trajectory_file"] = out.trajectory_file;
    d["final_file"] = out.final_file;
    d["trajectory_atoms"] = out.trajectory_atoms;
    d["lambdas"] = p.lambdas;

    py::list seq, dist, ang, wall, pos;
    for (const auto& r : p.restrseqs) seq.append(py::make_tuple(r.ai, r.aj, r.k, r.ih, r.to_center));
    for (const auto& r : p.restrdists) dist.append(py::make_tuple(r.ai, r.aj, r.d1, r.d2, r.k, r.ipsi));
    for (const auto& r : p.restrangs) ang.append(py::make_tuple(r.ai, r.aj, r.ak, r.ipsi, r.ang, r.k));
    for (const auto& r : p.restrwalls) wall.append(py::make_tuple(r.ai, r.aj, r.d, r.k, r.aMorse, r.dMorse, r.ih));
    for (const auto& r : p.restrspos)
        pos.append(py::make_tuple(r.a, r.ipsi, py::make_tuple(r.x.x, r.x.y, r.x.z), py::make_tuple(r.k.x, r.k.y, r.k.z)));
    d["seq_restraints"] = seq;
    d["dist_restraints"] = dist;
    d["ang_restraints"] = ang;
    d["wall_restraints"] = wall;
    d["pos_restraints"] = pos;
    return d;
}

// Config.from_inp: every setting of the .inp file. Paths come back absolute.
py::dict read_inp(const std::string& input_file) {
    InpParser parser(abs_inp(input_file));
    py::dict d = settings_dict(parser.parse_stage());
    d["topology"] = parser.topology_file();
    d["fep"] = parser.fep_file();
    d["restart"] = parser.restart_file();
    return d;
}

// Config(...): fill in what the user left unset, exactly as an .inp file with only
// these keys would. Derived defaults (charge_correction = polarisation, ...) included.
py::dict resolve_settings(const InpParser::Sections& sections, const InpParser::Lists& lists) {
    return settings_dict(InpParser::from_settings(sections, lists)->parse_stage());
}

class Engine {
   public:
    // Loads the system only (qgpu.Config's fixed settings). The first configure() sets the
    // run settings and, when initial_temperature is given, generates velocities from it and the seed.
    Engine(const std::string& topology, const std::string& fep,
           const InpParser::Sections& sections, const InpParser::Lists& lists, int device)
        : device_(device) {
        select_device(device_);
        const ParseResult parsed = InpParser::from_sections(topology, fep, sections, lists)->parse();
        handler_ = std::make_unique<CudaHandler>();
        handler_->initialize_system(parsed);

        base_stage_ = stage_from(parsed);  // template only: Python overwrites the run settings
        n_atoms_ = static_cast<int>(parsed.coords_init.size());
        n_states_ = static_cast<int>(parsed.lambdas.size());
    }

    StageSettings base_stage() const { return base_stage_; }
    int n_atoms() const { return n_atoms_; }
    int n_states() const { return n_states_; }

    void configure(const StageSettings& stage) {
        guard();
        check_cuda(cudaSetDevice(device_));
        handler_->validate_stage(stage);  // a rejected stage leaves the engine usable
        call([&] { handler_->configure(stage); });
    }

    void run() {
        guard();
        check_cuda(cudaSetDevice(device_));
        call([&] {
            handler_->run();
            check_cuda(cudaDeviceSynchronize());
        });
    }

    // Keeps the GIL: creates NumPy arrays.
    py::dict get_state() {
        guard();
        check_cuda(cudaSetDevice(device_));
        SystemState s;
        call([&] { s = handler_->get_state(); });

        const py::ssize_t n = static_cast<py::ssize_t>(s.coords.size());
        py::array_t<double> coords({n, py::ssize_t(3)});
        py::array_t<double> vels({n, py::ssize_t(3)});
        auto c = coords.mutable_unchecked<2>();
        auto v = vels.mutable_unchecked<2>();
        for (py::ssize_t i = 0; i < n; ++i) {
            c(i, 0) = s.coords[i].x;
            c(i, 1) = s.coords[i].y;
            c(i, 2) = s.coords[i].z;
            v(i, 0) = s.velocities[i].x;
            v(i, 1) = s.velocities[i].y;
            v(i, 2) = s.velocities[i].z;
        }
        py::dict out;
        out["coords"] = coords;
        out["velocities"] = vels;
        out["theta_corr"] = py::array_t<double>(s.theta_corr.size(), s.theta_corr.data());
        return out;
    }

    void set_state(py::array_t<double, py::array::c_style | py::array::forcecast> coords,
                   py::array_t<double, py::array::c_style | py::array::forcecast> velocities,
                   std::vector<double> theta_corr) {
        if (coords.ndim() != 2 || coords.shape(1) != 3 || coords.shape(0) != n_atoms_)
            throw InputError("coords must have shape (n_atoms, 3)");
        if (velocities.ndim() != 2 || velocities.shape(1) != 3 || velocities.shape(0) != n_atoms_)
            throw InputError("velocities must have shape (n_atoms, 3)");

        SystemState s;
        s.coords.resize(n_atoms_);
        s.velocities.resize(n_atoms_);
        auto c = coords.unchecked<2>();
        auto v = velocities.unchecked<2>();
        for (int i = 0; i < n_atoms_; ++i) {
            s.coords[i].x = c(i, 0);
            s.coords[i].y = c(i, 1);
            s.coords[i].z = c(i, 2);
            s.velocities[i].x = v(i, 0);
            s.velocities[i].y = v(i, 1);
            s.velocities[i].z = v(i, 2);
        }
        s.theta_corr = std::move(theta_corr);

        check_cuda(cudaSetDevice(device_));
        handler_->set_state(s);  // also clears the "broken" flag
        broken_ = false;
    }

   private:
    void guard() const {
        if (broken_)
            throw QgpuError("Engine is in a failed state after an error; call set_state() or create a new Engine.");
    }

    template <class F>
    void call(F&& f) {
        try {
            f();
        } catch (...) {
            broken_ = true;
            throw;
        }
    }

    int device_;
    int n_atoms_ = 0;
    int n_states_ = 0;
    bool broken_ = false;
    StageSettings base_stage_;
    std::unique_ptr<CudaHandler> handler_;
};

}  // namespace

// Built twice by CMake: qgpu._core (double) and qgpu._core_spfp (QDYN_SPFP).
#ifndef QGPU_MODULE_NAME
#define QGPU_MODULE_NAME _core
#endif

PYBIND11_MODULE(QGPU_MODULE_NAME, m) {
    m.doc() = "Private QGPU native module.";
#ifdef QDYN_SPFP
    m.attr("PRECISION") = "spfp";
#else
    m.attr("PRECISION") = "double";
#endif

    // Exceptions: raise the shared classes from qgpu/errors.py, so `except qgpu.QgpuError`
    // works for both precision modules. Leaked on purpose: no py::object destructors at exit.
    static const auto* errors = new py::module_(py::module_::import("qgpu.errors"));
    static const auto* qgpu_error = new py::object(errors->attr("QgpuError"));
    static const auto* input_error = new py::object(errors->attr("InputError"));
    static const auto* constraint_error = new py::object(errors->attr("ConstraintError"));
    static const auto* cuda_error = new py::object(errors->attr("CudaError"));
    py::register_exception_translator([](std::exception_ptr p) {
        try {
            if (p) std::rethrow_exception(p);
        } catch (const InputError& e) {
            py::set_error(*input_error, e.what());
        } catch (const ConstraintError& e) {
            py::set_error(*constraint_error, e.what());
        } catch (const CudaError& e) {
            py::set_error(*cuda_error, e.what());
        } catch (const QgpuError& e) {  // base class last
            py::set_error(*qgpu_error, e.what());
        }
    });

    // All classes are module_local: _core and _core_spfp register the same C++ types,
    // which pybind11 rejects as duplicates otherwise. pybind11 still converts these objects
    // between the two modules (same C++ type). That is only safe while every bound struct
    // has the same layout in both builds: never expose a type that contains real_t.

    // Plain structs
    py::class_<coord_t>(m, "Coord", py::module_local())
        .def(py::init<>())
        .def(py::init([](double x, double y, double z) { coord_t c; c.x = x; c.y = y; c.z = z; return c; }))
        .def_readwrite("x", &coord_t::x)
        .def_readwrite("y", &coord_t::y)
        .def_readwrite("z", &coord_t::z);

    py::class_<md_t>(m, "MD", py::module_local())  // stage fields only; System fields come from base_stage()
        .def_readwrite("steps", &md_t::steps)
        .def_readwrite("stepsize", &md_t::stepsize)
        .def_readwrite("temperature", &md_t::temperature)
        .def_readwrite("thermostat", &md_t::thermostat)
        .def_readwrite("bath_coupling", &md_t::bath_coupling)
        .def_readwrite("random_seed", &md_t::random_seed)
        .def_readwrite("initial_temperature", &md_t::initial_temperature)
        .def_readwrite("constraint_solvent", &md_t::constraint_solvent)
        .def_readwrite("constraint_solute", &md_t::constraint_solute)
        .def_readwrite("constraint_hydrogens", &md_t::constraint_hydrogens)
        .def_readwrite("shake_solvent", &md_t::shake_solvent)
        .def_readwrite("shake_solute", &md_t::shake_solute)
        .def_readwrite("shake_hydrogens", &md_t::shake_hydrogens)
        .def_readwrite("separate_scaling", &md_t::separate_scaling)
        .def_readwrite("shell_force", &md_t::shell_force)
        .def_readwrite("output", &md_t::output)
        .def_readwrite("energy", &md_t::energy)
        .def_readwrite("trajectory", &md_t::trajectory);

    py::class_<NativeOutputConfig>(m, "Outputs", py::module_local())
        .def_readwrite("final_file", &NativeOutputConfig::final_file)
        .def_readwrite("trajectory_file", &NativeOutputConfig::trajectory_file)
        .def_readwrite("energy_file", &NativeOutputConfig::energy_file)
        .def_readwrite("trajectory_atoms", &NativeOutputConfig::trajectory_atoms)
        .def_readonly("topology_file", &NativeOutputConfig::topology_file);

    // Restraints: atom indices are 1-based, exactly as in .inp files.
    py::class_<restrseq_t>(m, "SeqRestraint", py::module_local())
        .def(py::init([](int ai, int aj, double k, bool ih, int to_center) {
                 return restrseq_t{ai, aj, k, ih, to_center};
             }),
             py::arg("ai"), py::arg("aj"), py::arg("k"), py::arg("ih") = false, py::arg("to_center") = 0);
    py::class_<restrpos_t>(m, "PosRestraint", py::module_local())
        .def(py::init([](int a, int ipsi, coord_t x, coord_t k) { return restrpos_t{a, ipsi, x, k}; }),
             py::arg("a"), py::arg("ipsi"), py::arg("x"), py::arg("k"));
    py::class_<restrdis_t>(m, "DistRestraint", py::module_local())
        .def(py::init([](int ai, int aj, double d1, double d2, double k, int ipsi) {
                 restrdis_t r{};
                 r.ai = ai;
                 r.aj = aj;
                 r.d1 = d1;
                 r.d2 = d2;
                 r.k = k;
                 r.ipsi = ipsi;
                 return r;
             }),
             py::arg("ai"), py::arg("aj"), py::arg("d1"), py::arg("d2"), py::arg("k"), py::arg("ipsi") = 0);
    py::class_<restrang_t>(m, "AngRestraint", py::module_local())
        .def(py::init([](int ai, int aj, int ak, int ipsi, double ang, double k) {
                 return restrang_t{ai, aj, ak, ipsi, ang, k};
             }),
             py::arg("ai"), py::arg("aj"), py::arg("ak"), py::arg("ipsi"), py::arg("ang"), py::arg("k"));
    py::class_<restrwall_t>(m, "WallRestraint", py::module_local())
        .def(py::init([](int ai, int aj, double d, double k, double a_morse, double d_morse, bool ih) {
                 return restrwall_t{ai, aj, d, k, a_morse, d_morse, ih};
             }),
             py::arg("ai"), py::arg("aj"), py::arg("d"), py::arg("k"),
             py::arg("a_morse"), py::arg("d_morse"), py::arg("ih") = false);

    py::class_<StageSettings>(m, "StageSettings", py::module_local())
        .def(py::init<const StageSettings&>())  // copy
        .def_readwrite("md", &StageSettings::md)
        .def_readwrite("lambdas", &StageSettings::lambdas)
        .def_readwrite("restrspos", &StageSettings::restrspos)
        .def_readwrite("restrseqs", &StageSettings::restrseqs)
        .def_readwrite("restrdists", &StageSettings::restrdists)
        .def_readwrite("restrangs", &StageSettings::restrangs)
        .def_readwrite("restrwalls", &StageSettings::restrwalls)
        .def_readwrite("outputs", &StageSettings::outputs);

    py::class_<Engine>(m, "Engine", py::module_local())
        .def(py::init<const std::string&, const std::string&, const InpParser::Sections&,
                      const InpParser::Lists&, int>(),
             py::arg("topology"), py::arg("fep"), py::arg("sections"), py::arg("lists"), py::arg("device") = 0)
        .def_property_readonly("n_atoms", &Engine::n_atoms)
        .def_property_readonly("n_states", &Engine::n_states)
        .def("base_stage", &Engine::base_stage)
        .def("configure", &Engine::configure, py::arg("stage"), py::call_guard<py::gil_scoped_release>())
        .def("run", &Engine::run, py::call_guard<py::gil_scoped_release>())
        .def("get_state", &Engine::get_state)
        .def("set_state", &Engine::set_state, py::arg("coords"), py::arg("velocities"), py::arg("theta_corr"));

    m.def("read_inp", &read_inp, py::arg("input_file"));
    m.def("resolve_settings", &resolve_settings, py::arg("sections"), py::arg("lists"));
}