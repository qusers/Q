"""Setup script to build and install the QligFEP package with Fortran compilation.
It integrates the Fortran compilation step into the setuptools build process.
"""

import platform
import shutil
import subprocess
import sys
from pathlib import Path

from setuptools import setup
from setuptools.command.build_py import build_py as _build_py


def find_gfortran():
    """Find gfortran compiler in environment or system PATH."""
    gfortran = shutil.which("gfortran")
    if gfortran:
        return gfortran

    if platform.system() == "Darwin":
        homebrew_paths = [
            "/opt/homebrew/bin/gfortran",
            "/usr/local/bin/gfortran",
        ]
        for path in homebrew_paths:
            if Path(path).exists():
                return path

    return None


def compile_qprep():
    """Compile qprep binary using make."""
    src_dir = Path(__file__).parent / "src"
    q6_dir = src_dir / "q6"

    if not q6_dir.exists():
        print(f"WARNING: Q6 source directory not found: {q6_dir}", file=sys.stderr)
        return False

    gfortran = find_gfortran()
    if not gfortran:
        system = platform.system()
        comp_flag = "osx" if system == "Darwin" else "gcc"
        print(
            f"\nWARNING: gfortran compiler not found!\n"
            f"The qprep binary could not be compiled automatically.\n"
            f"Please install gfortran and manually compile:\n"
            f"  cd {q6_dir}\n"
            f"  make qprep COMP={comp_flag}\n",
            file=sys.stderr
        )
        return False

    print(f"\nFound gfortran: {gfortran}")
    print(f"Compiling qprep in {q6_dir}...")

    system = platform.system()
    comp_flag = "osx" if system == "Darwin" else "gcc"
    print(f"Using compiler configuration: COMP={comp_flag}\n")

    try:
        subprocess.run(
            ["make", "clean"],
            cwd=q6_dir,
            check=False,
            capture_output=True
        )

        result = subprocess.run(
            ["make", "qprep", f"COMP={comp_flag}"],
            cwd=q6_dir,
            check=True,
            capture_output=True,
            text=True
        )

        print(result.stdout)

        bin_dir = q6_dir / "bin" / "q6"
        bin_dir.mkdir(parents=True, exist_ok=True)

        qprep_binary = q6_dir / "qprep"
        if qprep_binary.exists():
            shutil.move(str(qprep_binary), str(bin_dir / "qprep"))
            print(f"Successfully compiled qprep -> {bin_dir / 'qprep'}\n")

            obj_dir = q6_dir / "obj"
            obj_dir.mkdir(exist_ok=True)
            for pattern in ["*.o", "*.mod"]:
                for obj_file in q6_dir.glob(pattern):
                    shutil.move(str(obj_file), str(obj_dir / obj_file.name))

            return True
        else:
            print("ERROR: qprep binary not found after compilation", file=sys.stderr)
            return False

    except subprocess.CalledProcessError as e:
        print(
            f"\nERROR: Compilation failed!\n"
            f"Command: {' '.join(e.cmd)}\n"
            f"Output: {e.stdout}\n"
            f"Error: {e.stderr}",
            file=sys.stderr
        )
        return False


class BuildWithFortran(_build_py):
    """Custom build_py command that compiles Fortran code before building Python package."""

    def run(self):
        """Execute the build, including Fortran compilation."""
        print("\n" + "="*70)
        print("Building QligFEP - Compiling Fortran qprep binary")
        print("="*70)

        compile_qprep()

        print("="*70)
        print("Continuing with Python package build")
        print("="*70 + "\n")

        super().run()


setup(
    cmdclass={
        'build_py': BuildWithFortran,
    }
)
