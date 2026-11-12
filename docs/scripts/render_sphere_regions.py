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

import numpy as np
from PIL import Image, ImageFilter
from pymol import cgo, cmd

OUT = Path(__file__).resolve().parents[1] / "source" / "_static" / "figures" / "sphere-regions.png"
NORMAL = (0.0, 1.0, 0.0)  # the slab is perpendicular to this axis; the camera looks along it
UP = (0.0, 0.0, 1.0)
HALF_SLAB = 1.8  # A on each side of the plane through the centre
STICK_RADIUS = 0.13  # A, protein, water and ligand
FREE = (0.10, 0.55, 0.70)  # blue
EXCLUDED = (0.79, 0.76, 0.70)  # pale sand; also the colour of its surface
WATER_O = (0.88, 0.12, 0.10)
WATER_H = (0.94, 0.95, 0.96)
LIGAND = (0.98, 0.50, 0.08)
CIRCLE = (0.46, 0.48, 0.51)  # mid grey, readable on white and on the dark theme
EDGE = (69, 74, 80)  # one-pixel silhouette that keeps the white hydrogens visible on white
SIZE = 1500  # px


def dashed_circle(radius, dashes=64, width=0.085):
    """Dashed circle of the sphere radius in the plane z = 0 of the rotated frame."""
    obj = []
    for i in range(dashes):
        start = 2.0 * math.pi * i / dashes
        end = start + 0.56 * 2.0 * math.pi / dashes
        angles = np.linspace(start, end, 5)
        for a, b in zip(angles[:-1], angles[1:]):
            p = (radius * math.cos(a), radius * math.sin(a), 0.0)
            q = (radius * math.cos(b), radius * math.sin(b), 0.0)
            obj += [cgo.CYLINDER, *p, *q, width, *CIRCLE, *CIRCLE]
    return obj


def ray(path: Path) -> Image.Image:
    cmd.png(str(path), width=SIZE, height=SIZE, dpi=200, ray=1)
    return Image.open(path).convert("RGBA")


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--benchmark-dir", type=Path, required=True)
    parser.add_argument("--free", type=Path, required=True)
    args = parser.parse_args()
    prep = json.loads((args.benchmark_dir / "preparation.json").read_text())
    centre, radius = np.array(prep["legacy_cog"]), float(prep["nominal_radius"])

    # Rotate everything into a frame with the sphere centre at the origin and the slab normal along z.
    normal = np.array(NORMAL) / np.linalg.norm(NORMAL)
    u = np.cross(UP, normal)
    u /= np.linalg.norm(u)
    rotation = np.array([u, np.cross(normal, u), normal])
    inputs = {
        "solute": args.benchmark_dir / "protein.pdb",
        "mask": args.free,
        "water": args.benchmark_dir / "water.pdb",
        "ligand": args.benchmark_dir / "ejm_31.pdb",
    }
    for name, path in inputs.items():
        cmd.load(str(path), name)
    cmd.remove("hydrogens and not water")
    # A solute residue is free if qprep kept it in the mask "not excluded".
    cmd.select("free_solute", "solute and (mask around 0.01)")
    cmd.select("excluded_solute", "solute and not free_solute")
    for name in ("free_solute", "excluded_solute"):
        print(f"{name}: {cmd.count_atoms(name)} heavy atoms, {cmd.count_atoms(f'{name} and name CA')} residues")
    for name in inputs:
        cmd.load_coords((cmd.get_coords(name) - centre) @ rotation.T, name)

    slab = f"z > {-HALF_SLAB} and z < {HALF_SLAB}"
    cmd.select("solute_slice", f"byres (solute and {slab})")
    cmd.select("water_slice", f"byres (water and elem O and {slab})")
    print(f"drawn: {cmd.count_atoms('solute_slice and free_solute and name CA')} free and "
          f"{cmd.count_atoms('solute_slice and excluded_solute and name CA')} excluded residues, "
          f"{cmd.count_atoms('water_slice and elem O')} water molecules")

    cmd.hide("everything")
    colours = {"c_free": FREE, "c_excluded": EXCLUDED, "c_water_o": WATER_O, "c_water_h": WATER_H, "c_ligand": LIGAND}
    for name, rgb in colours.items():
        cmd.set_color(name, list(rgb))
    cmd.color("c_free", "free_solute")
    cmd.color("c_excluded", "excluded_solute")
    cmd.color("c_water_o", "water and elem O")
    cmd.color("c_water_h", "water and elem H")
    cmd.color("c_ligand", "ligand")
    cmd.show("sticks", "solute_slice or water_slice or ligand")
    cmd.set("stick_radius", STICK_RADIUS)
    cmd.set("stick_ball", 0)
    cmd.set("valence", 0)
    # Translucent surface around the excluded residues that are drawn.
    cmd.create("excluded_shell", "excluded_solute and solute_slice")
    cmd.hide("everything", "excluded_shell")
    cmd.show("surface", "excluded_shell")
    cmd.set("transparency", 0.78, "excluded_shell")

    cmd.set("orthoscopic", 1)
    cmd.set("depth_cue", 0)
    cmd.set("ray_shadows", 0)
    cmd.set("specular", 0.08)
    cmd.set("ambient", 0.65)
    cmd.set("direct", 0.35)
    cmd.set("ray_opaque_background", 0)
    cmd.set("antialias", 2)
    # Frame the drawn atoms and the whole circle.
    xy = cmd.get_coords("solute_slice or ligand or water_slice")[:, :2]
    low = np.minimum(xy.min(0) - 2.8, -radius - 1.4)
    high = np.maximum(xy.max(0) + 2.8, radius + 1.4)
    middle = (low + high) / 2.0
    cmd.set_view([1, 0, 0, 0, 1, 0, 0, 0, 1, 0, 0, -200, *middle, 0, 100, 300, 0])
    cmd.pseudoatom("frame", pos=[*middle, 0.0])
    cmd.zoom("frame", float(max(high - low) / 2.0))
    cmd.hide("everything", "frame")
    view = cmd.get_view()

    scratch = Path(tempfile.gettempdir())
    molecules = ray(scratch / "sphere-regions-molecules.png")
    edge = Image.new("RGBA", molecules.size, (*EDGE, 0))
    edge.putalpha(molecules.getchannel("A").filter(ImageFilter.MaxFilter(3)))
    molecules = Image.alpha_composite(edge, molecules)

    # The circle is rendered on its own in the same view, so the molecules cannot hide it.
    cmd.hide("everything")
    cmd.load_cgo(dashed_circle(radius), "circle")
    cmd.set_view(view)
    circle = ray(scratch / "sphere-regions-circle.png")

    image = Image.alpha_composite(molecules, circle)
    image = image.crop(image.getbbox())
    # 256 colours with alpha keep the file small; the render has few colours
    image.quantize(colors=256, method=Image.Quantize.FASTOCTREE, dither=Image.Dither.NONE).save(OUT, optimize=True)
    print(f"wrote {OUT} ({image.width} x {image.height} px, {OUT.stat().st_size / 1e6:.2f} MB)")


if __name__ == "__main__":
    main()
