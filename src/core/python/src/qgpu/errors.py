"""QGPU exceptions, shared by every native precision module (_core, _core_spfp)."""


class QgpuError(RuntimeError):
    """Base class for errors raised by the native QGPU engine."""


class InputError(QgpuError):
    """Invalid input: files, settings, stages or state."""


class ConstraintError(QgpuError):
    """A constraint solver (SHAKE/LINCS/SETTLE) failed."""


class CudaError(QgpuError):
    """A CUDA call failed."""
