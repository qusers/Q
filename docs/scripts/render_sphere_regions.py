"""Figure for docs/source/concepts/spherical-boundary.md: regions of the simulation sphere.

Renders a slab through the centre of the Tyk2 sphere of the benchmark set with
PyMOL. The structures are read in place and are not part of the repository.

Inputs
- ``--benchmark-dir``: the ``tyk2`` folder with ``protein.pdb``, ``water.pdb``,
  ``ejm_31.pdb`` and ``preparation.json`` (sphere centre and radius).
- ``--free``: PDB of the solute atoms that are not excluded, written by qprep:

      rl AMBER14sb.lib
      rprm AMBER14sb.prm
      rp protein.pdb
      boundary 1 <x> <y> <z> 25.0
      maketop tyk2
      mask none
      mask not excluded
      wp notexcluded.pdb y

    python docs/scripts/render_sphere_regions.py \
        --benchmark-dir /path/to/perturbations/tyk2 --free /path/to/notexcluded.pdb
"""

from __future__ import annotations

import argparse
import json
import math
import tempfile
from pathlib import Path

from PIL import Image
from pymol import cgo, cmd

OUT = Path(__file__).resolve().parents[1] / "source" / "_static" / "figures" / "sphere-regions.png"
HALF_SLAB = 3.0  # A on each side of the plane through the centre
FREE = (0.10, 0.62, 0.56)  # teal
EXCLUDED = (0.60, 0.60, 0.60)  # grey
WATER = (0.22, 0.53, 0.90)  # blue
LIGAND = (0.92, 0.41, 0.20)  # orange
OUTLINE = (0.50, 0.50, 0.50)  # readable on white and on the dark theme
SIZE = 1500  # px


# Slab normal and the two in-plane axes that become the image x and y directions
NORMAL, U, V = 0, 1, 2  # slab perpendicular to x; image axes y and z
AXES = "xyz"


def point(centre, a, b):
    """Point in the slab plane at in-plane offsets a (image x) and b (image y)."""
    p = list(centre)
    p[U] += a
    p[V] += b
    return p


def dashed_circle(centre, radius, dashes=72, width=0.18):
    """Dashed circle in the slab plane through the centre, as CGO cylinders."""
    obj = []
    for i in range(dashes):
        a0 = 2.0 * math.pi * i / dashes
        a1 = a0 + 0.55 * 2.0 * math.pi / dashes
        p0 = point(centre, radius * math.cos(a0), radius * math.sin(a0))
        p1 = point(centre, radius * math.cos(a1), radius * math.sin(a1))
        obj += [cgo.CYLINDER, *p0, *p1, width, *OUTLINE, *OUTLINE]
    return obj


def cross(centre, arm=1.2, width=0.18):
    return [
        cgo.CYLINDER, *point(centre, -arm, 0), *point(centre, arm, 0), width, *OUTLINE, *OUTLINE,
        cgo.CYLINDER, *point(centre, 0, -arm), *point(centre, 0, arm), width, *OUTLINE, *OUTLINE,
    ]


def view_rotation():
    """PyMOL view rotation (column-major) that maps U, V and the normal to the screen x, y and z."""
    rotation = [0.0] * 9
    for screen, model in enumerate((U, V, NORMAL)):
        rotation[model * 3 + screen] = 1.0
    return rotation


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--benchmark-dir", type=Path, required=True)
    parser.add_argument("--free", type=Path, required=True)
    args = parser.parse_args()
    prep = json.loads((args.benchmark_dir / "preparation.json").read_text())
    centre, radius = prep["legacy_cog"], prep["nominal_radius"]

    cmd.load(str(args.benchmark_dir / "protein.pdb"), "protein")
    cmd.load(str(args.free), "free")
    cmd.load(str(args.benchmark_dir / "water.pdb"), "water")
    cmd.load(str(args.benchmark_dir / "ejm_31.pdb"), "ligand")
    cmd.remove("hydrogens")

    n = AXES[NORMAL]
    slab = f"({n} > {centre[NORMAL] - HALF_SLAB} and {n} < {centre[NORMAL] + HALF_SLAB})"
    # A solute residue is free if qprep kept it in the mask "not excluded".
    cmd.select("free_solute", "protein and (free around 0.01)")
    cmd.select("excluded_solute", "protein and not free_solute")
    for name in ("free_solute", "excluded_solute"):
        print(f"{name}: {cmd.count_atoms(name)} heavy atoms, {cmd.count_atoms(f'{name} and name CA')} residues")

    cmd.hide("everything")
    cmd.show("sticks", f"byres (protein and {slab})")
    cmd.show("spheres", f"water and {slab}")
    cmd.show("sticks", "ligand")
    cmd.set("stick_radius", 0.22)
    cmd.set("stick_radius", 0.45, "ligand")
    cmd.set("sphere_scale", 0.32, "water")
    cmd.set_color("c_free", FREE)
    cmd.set_color("c_excluded", EXCLUDED)
    cmd.set_color("c_water", WATER)
    cmd.set_color("c_ligand", LIGAND)
    cmd.color("c_free", "free_solute")
    cmd.color("c_excluded", "excluded_solute")
    cmd.color("c_water", "water")
    cmd.color("c_ligand", "ligand")
    cmd.load_cgo(dashed_circle(centre, radius) + cross(centre), "outline")

    # Look along the slab normal at the centre, without perspective, depth cue or shadows.
    cmd.set("orthoscopic", 1)
    cmd.set("depth_cue", 0)
    cmd.set("ray_shadows", 0)
    cmd.set("specular", 0)
    cmd.set("ambient", 0.6)
    cmd.set("ray_opaque_background", 0)
    cmd.set("antialias", 2)
    cmd.set_view([*view_rotation(), 0, 0, -200, *centre, 150, 250, 0])
    cmd.pseudoatom("frame", pos=centre)
    cmd.zoom("frame", 42)
    cmd.hide("everything", "frame")
    raw = Path(tempfile.gettempdir()) / "sphere-regions-raw.png"
    cmd.png(str(raw), width=SIZE, height=SIZE, dpi=200, ray=1)

    image = Image.open(raw).convert("RGBA")
    image = image.crop(image.getbbox())
    # 256 colours with alpha keep the file small; the render has few colours
    image.quantize(colors=256, method=Image.Quantize.FASTOCTREE, dither=Image.Dither.NONE).save(OUT, optimize=True)
    print(f"wrote {OUT} ({image.width} x {image.height} px, {OUT.stat().st_size / 1e6:.2f} MB)")


if __name__ == "__main__":
    main()
