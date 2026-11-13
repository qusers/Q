"""Restraint types. Plain Python, so they work with every precision module.

Atom indices are 1-based, exactly as in .inp files. Each class converts itself
to the matching native struct of the engine's module with `_native(module)`.
"""

from dataclasses import dataclass


@dataclass(frozen=True)
class Coord:
    x: float = 0.0
    y: float = 0.0
    z: float = 0.0

    def _native(self, mod):
        return mod.Coord(self.x, self.y, self.z)


@dataclass(frozen=True)
class SeqRestraint:
    ai: int
    aj: int
    k: float
    ih: bool = False
    to_center: int = 0

    def _native(self, mod):
        return mod.SeqRestraint(self.ai, self.aj, self.k, self.ih, self.to_center)


@dataclass(frozen=True)
class PosRestraint:
    a: int
    ipsi: int
    x: Coord
    k: Coord

    def __post_init__(self):
        object.__setattr__(self, "x", _coord(self.x))
        object.__setattr__(self, "k", _coord(self.k))

    def _native(self, mod):
        return mod.PosRestraint(self.a, self.ipsi, self.x._native(mod), self.k._native(mod))


@dataclass(frozen=True)
class DistRestraint:
    ai: int
    aj: int
    d1: float
    d2: float
    k: float
    ipsi: int = 0

    def _native(self, mod):
        return mod.DistRestraint(self.ai, self.aj, self.d1, self.d2, self.k, self.ipsi)


@dataclass(frozen=True)
class AngRestraint:
    ai: int
    aj: int
    ak: int
    ipsi: int
    ang: float
    k: float

    def _native(self, mod):
        return mod.AngRestraint(self.ai, self.aj, self.ak, self.ipsi, self.ang, self.k)


@dataclass(frozen=True)
class WallRestraint:
    ai: int
    aj: int
    d: float
    k: float
    a_morse: float
    d_morse: float
    ih: bool = False

    def _native(self, mod):
        return mod.WallRestraint(self.ai, self.aj, self.d, self.k, self.a_morse, self.d_morse, self.ih)


def _coord(value):
    """Accept a Coord or any (x, y, z) sequence."""
    return value if isinstance(value, Coord) else Coord(*value)
