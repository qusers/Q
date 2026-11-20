"""Select the native extension for a floating-point precision."""

from importlib import import_module
from importlib.util import find_spec

# Precision -> extension module, built by setup.py (see QGPU_PRECISIONS).
MODULES = {"double": "_core", "spfp": "_core_spfp"}


def module_name(precision):
    if precision not in MODULES:
        raise ValueError(f"precision must be one of {sorted(MODULES)}, got {precision!r}")
    return f"qgpu.{MODULES[precision]}"


def is_available(precision):
    """True if the extension for `precision` is installed, without importing it."""
    return find_spec(module_name(precision)) is not None


def load(precision):
    name = module_name(precision)
    try:
        return import_module(name)
    except ModuleNotFoundError as e:
        if e.name != name:
            raise
        raise RuntimeError(
            f"qgpu was built without {precision!r} support. Reinstall with "
            f"QGPU_PRECISIONS including {precision!r}.") from e


def load_any():
    """Any installed module, for helpers that don't depend on precision (read_inp)."""
    for precision in MODULES:
        if is_available(precision):
            return load(precision)
    raise RuntimeError("No qgpu native extension is installed.")
