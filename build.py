"""Custom build script to compile Fortran qprep binary during pip installation.
It detects the gfortran compiler and runs make to build the qprep executable.
"""

import platform
import shutil
import subprocess
import sys
from pathlib import Path

from setuptools.command.build_py import build_py as _build_py


def find_gfortran():
    """
    Find gfortran compiler in environment or system PATH.

    Returns:
        Path to gfortran executable or None if not found
    """
    # Check if gfortran is available
    gfortran = shutil.which("gfortran")
    if gfortran:
        return gfortran

    # On macOS, check common homebrew locations
    if platform.system() == "Darwin":
        homebrew_paths = [
            "/opt/homebrew/bin/gfortran",  # Apple Silicon
            "/usr/local/bin/gfortran",      # Intel Mac
        ]
        for path in homebrew_paths:
            if Path(path).exists():
                return path

    return None


def compile_qprep(src_dir):
    """
    Compile qprep binary using make.

    Args:
        src_dir: Path to the Q source directory
    """
    q6_dir = src_dir / "q6"

    if not q6_dir.exists():
        raise FileNotFoundError(f"Q6 source directory not found: {q6_dir}")

    # Check for gfortran
    gfortran = find_gfortran()
    if not gfortran:
        # Determine the appropriate COMP flag for the platform
        system = platform.system()
        comp_flag = "osx" if system == "Darwin" else "gcc"

        print(
            "\n" + "="*70,
            "\nWARNING: gfortran compiler not found!",
            "\n" + "="*70,
            "\nThe qprep binary could not be compiled automatically.",
            "\nYou have two options:",
            "\n  1. Install gfortran and reinstall this package:",
            "\n     - On macOS: micromamba install gfortran -c conda-forge",
            "\n     - On Linux: micromamba install gfortran -c conda-forge",
            "\n  2. Manually compile qprep after installation:",
            f"\n     cd {q6_dir}",
            f"\n     make qprep COMP={comp_flag}",
            "\n" + "="*70 + "\n",
            file=sys.stderr
        )
        return False

    print(f"Found gfortran: {gfortran}")
    print(f"Compiling qprep in {q6_dir}...")

    # Determine the appropriate COMP flag for the platform
    system = platform.system()
    comp_flag = "osx" if system == "Darwin" else "gcc"
    print(f"Using compiler configuration: COMP={comp_flag}")

    try:
        # Clean any previous build artifacts
        subprocess.run(
            ["make", "clean"],
            cwd=q6_dir,
            check=False,  # Don't fail if clean fails
            capture_output=True
        )

        # Compile qprep only
        result = subprocess.run(
            ["make", "qprep", f"COMP={comp_flag}"],
            cwd=q6_dir,
            check=True,
            capture_output=True,
            text=True
        )

        print("Compilation output:")
        print(result.stdout)

        # Create bin directory structure
        bin_dir = q6_dir / "bin" / "q6"
        bin_dir.mkdir(parents=True, exist_ok=True)

        # Move qprep binary to bin directory
        qprep_binary = q6_dir / "qprep"
        if qprep_binary.exists():
            shutil.move(str(qprep_binary), str(bin_dir / "qprep"))
            print(f"Successfully compiled qprep -> {bin_dir / 'qprep'}")

            # Move object files to obj directory
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
            f"Return code: {e.returncode}\n"
            f"Output: {e.stdout}\n"
            f"Error: {e.stderr}",
            file=sys.stderr
        )
        return False


def build():
    """Main build function called by setuptools."""
    # Get the source directory
    src_dir = Path(__file__).parent / "src"

    print("\n" + "="*70)
    print("Building QligFEP - Compiling Fortran qprep binary")
    print("="*70 + "\n")

    success = compile_qprep(src_dir)

    if not success:
        print(
            "\nWARNING: Package installed but qprep compilation failed.",
            "\nYou will need to compile qprep manually before using QligFEP.",
            file=sys.stderr
        )
        # Don't fail the installation - let users install and compile manually
        # This is important for HPC environments where they may need specific compilers
    else:
        print("\n" + "="*70)
        print("QligFEP build completed successfully!")
        print("="*70 + "\n")


class BuildWithFortran(_build_py):
    """Custom build_py command that compiles Fortran code before building Python package."""

    def run(self):
        """Execute the build, including Fortran compilation."""
        # Compile Fortran first
        build()
        # Then run the standard build_py
        super().run()


if __name__ == "__main__":
    build()
