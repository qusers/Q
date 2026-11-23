import os
from pathlib import Path
import subprocess
import sys

from setuptools import Extension, setup
from setuptools.command.build_ext import build_ext


PACKAGE_ROOT = Path(__file__).resolve().parent
CORE_ROOT = PACKAGE_ROOT.parent

# Precision -> extension module (and CMake target) name.
MODULES = {"double": "_core", "spfp": "_core_spfp"}


class CMakeExtension(Extension):
    def __init__(self, name, target):
        super().__init__(name, sources=[])
        self.target = target


class CMakeBuild(build_ext):
    def build_extension(self, ext):
        import pybind11

        if not (CORE_ROOT / "CMakeLists.txt").is_file():
            raise RuntimeError(
                "QGPU must be built from the full repository checkout: "
                f"missing {CORE_ROOT / 'CMakeLists.txt'}"
            )

        extension_path = Path(self.get_ext_fullpath(ext.name)).resolve()
        output_dir = extension_path.parent
        build_dir = Path(self.build_temp).resolve() / "cmake"  # shared by all extensions

        output_dir.mkdir(parents=True, exist_ok=True)
        build_dir.mkdir(parents=True, exist_ok=True)

        config = "Debug" if self.debug else "Release"
        architecture = os.environ.get("QGPU_CUDA_ARCH", "89")
        parallel = os.environ.get("CMAKE_BUILD_PARALLEL_LEVEL", "4")

        configure_command = [
            "cmake",
            "-S", str(CORE_ROOT),
            "-B", str(build_dir),
            f"-DCMAKE_BUILD_TYPE={config}",
            "-DQGPU_BUILD_PYTHON=ON",
            "-DQDYN_BUILD_TESTS=OFF",
            f"-DPython_EXECUTABLE={sys.executable}",
            f"-Dpybind11_DIR={pybind11.get_cmake_dir()}",
            f"-DCMAKE_CUDA_ARCHITECTURES={architecture}",
            "-DQDYN_SPFP=OFF",  # _core reuses the double core; _core_spfp gets its own
            f"-DCMAKE_LIBRARY_OUTPUT_DIRECTORY={output_dir}",
            f"-DCMAKE_LIBRARY_OUTPUT_DIRECTORY_{config.upper()}={output_dir}",
        ]

        subprocess.run(configure_command, check=True)

        subprocess.run(
            [
                "cmake",
                "--build", str(build_dir),
                "--config", config,
                "--target", ext.target,
                "--parallel", parallel,
            ],
            check=True,
        )

        if not extension_path.is_file():
            raise RuntimeError(
                "CMake completed, but the expected Python extension "
                f"was not produced: {extension_path}"
            )


def selected_precisions():
    names = [p.strip() for p in os.environ.get("QGPU_PRECISIONS", "double,spfp").split(",") if p.strip()]
    unknown = sorted(set(names) - MODULES.keys())
    if unknown or not names:
        raise RuntimeError(
            f"QGPU_PRECISIONS must list one or more of {sorted(MODULES)}, got {names or 'nothing'}")
    return list(dict.fromkeys(names))


setup(
    ext_modules=[CMakeExtension(f"qgpu.{MODULES[p]}", MODULES[p]) for p in selected_precisions()],
    cmdclass={"build_ext": CMakeBuild},
)
