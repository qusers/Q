"""Figure for docs/source/concepts/electrostatics.md: ionisable residues by treatment.

Renders the Tyk2 structure of the benchmark set with PyMOL, as prepared by
``qprep_prot``. The residues that were neutralised are read from
``preparation.json``; the structures are read in place and are not part of
the repository.

    python docs/scripts/render_boundary_neutralisation.py --benchmark-dir /path/to/perturbations/tyk2
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

OUT = Path(__file__).resolve().parents[1] / "source" / "_static" / "figures" / "boundary-neutralisation.png"
NORMAL = (1.0, 0.0, 0.0)  # same view as the sphere-regions figure
UP = (0.0, 0.0, 1.0)
HALF_SLAB = 8.0  # A; ionisable residues with their C-alpha in this slab are drawn
STICK_RADIUS = 0.2
ACIDIC = (0.85, 0.20, 0.20)  # ASP, GLU, charged
BASIC = (0.16, 0.45, 0.85)  # LYS, ARG, HIP, charged
NEUTRALISED = (0.79, 0.76, 0.70)  # ASH, GLH, LYN, ARN
BACKBONE = (0.55, 0.57, 0.60)
LIGAND = (0.98, 0.50, 0.08)
CIRCLE = (0.46, 0.48, 0.51)
EDGE = (69, 74, 80)
SIZE = 1500  # px
CHARGED_ACIDIC, CHARGED_BASIC = "ASP+GLU", "LYS+ARG+HIP"


def dashed_circle(radius, dashes=64, width=0.085):
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
    args = parser.parse_args()
    prep = json.loads((args.benchmark_dir / "preparation.json").read_text())
    centre, radius = np.array(prep["legacy_cog"]), float(prep["nominal_radius"])
    changed = prep["neutralization"]["modifications"]
    neutralised_ids = "+".join(key.split(":")[1] for key in changed)

    normal = np.array(NORMAL) / np.linalg.norm(NORMAL)
    u = np.cross(UP, normal)
    u /= np.linalg.norm(u)
    rotation = np.array([u, np.cross(normal, u), normal])
    cmd.load(str(args.benchmark_dir / "protein.pdb"), "solute")
    cmd.load(str(args.benchmark_dir / "ejm_31.pdb"), "ligand")
    cmd.remove("hydrogens")
    for name in ("solute", "ligand"):
        cmd.load_coords((cmd.get_coords(name) - centre) @ rotation.T, name)

    cmd.select("neutralised", f"solute and resi {neutralised_ids}")
    cmd.select("acidic", f"solute and resn {CHARGED_ACIDIC}")
    cmd.select("basic", f"solute and resn {CHARGED_BASIC}")
    in_slab = f"byres (solute and name CA and z > {-HALF_SLAB} and z < {HALF_SLAB})"
    for name in ("neutralised", "acidic", "basic"):
        distances = [np.linalg.norm(cmd.get_coords(f"solute and resi {a.resi}").mean(0))
                     for a in cmd.get_model(f"{name} and name CA").atom]
        beyond_22 = sum(d > radius - 3.0 for d in distances)
        print(f"{name}: {len(distances)} residues, {cmd.count_atoms(f'{name} and name CA and {in_slab}')} drawn, "
              f"centre distances {min(distances):.1f} to {max(distances):.1f} A, {beyond_22} beyond {radius - 3.0:.0f} A")

    cmd.hide("everything")
    colours = {"c_acidic": ACIDIC, "c_basic": BASIC, "c_neutral": NEUTRALISED, "c_backbone": BACKBONE, "c_ligand": LIGAND}
    for name, rgb in colours.items():
        cmd.set_color(name, list(rgb))
    cmd.color("c_backbone", "solute")
    cmd.color("c_acidic", "acidic")
    cmd.color("c_basic", "basic")
    cmd.color("c_neutral", "neutralised")
    cmd.color("c_ligand", "ligand")
    # Backbone trace of the residues in the slab, for context
    cmd.set("ribbon_width", 2.0)
    cmd.set("ribbon_sampling", 8)
    cmd.show("ribbon", f"solute and {in_slab}")
    cmd.show("sticks", f"(acidic or basic or neutralised) and {in_slab} and (sidechain or name CA)")
    cmd.show("sticks", "ligand")
    cmd.set("stick_radius", STICK_RADIUS)
    cmd.set("stick_radius", 0.13, "ligand")
    cmd.set("stick_ball", 0)
    cmd.set("valence", 0)

    cmd.set("orthoscopic", 1)
    cmd.set("depth_cue", 0)
    cmd.set("ray_shadows", 0)
    cmd.set("specular", 0.08)
    cmd.set("ambient", 0.65)
    cmd.set("direct", 0.35)
    cmd.set("ray_opaque_background", 0)
    cmd.set("antialias", 2)
    xy = cmd.get_coords(f"(solute and {in_slab}) or ligand")[:, :2]
    low = np.minimum(xy.min(0) - 2.8, -radius - 1.4)
    high = np.maximum(xy.max(0) + 2.8, radius + 1.4)
    middle = (low + high) / 2.0
    cmd.set_view([1, 0, 0, 0, 1, 0, 0, 0, 1, 0, 0, -200, *middle, 0, 100, 300, 0])
    cmd.pseudoatom("frame", pos=[*middle, 0.0])
    cmd.zoom("frame", float(max(high - low) / 2.0))
    cmd.hide("everything", "frame")
    view = cmd.get_view()

    scratch = Path(tempfile.gettempdir())
    molecules = ray(scratch / "boundary-neutralisation-molecules.png")
    edge = Image.new("RGBA", molecules.size, (*EDGE, 0))
    edge.putalpha(molecules.getchannel("A").filter(ImageFilter.MaxFilter(3)))
    molecules = Image.alpha_composite(edge, molecules)

    cmd.hide("everything")
    cmd.load_cgo(dashed_circle(radius), "circle")
    cmd.set_view(view)
    circle = ray(scratch / "boundary-neutralisation-circle.png")

    image = Image.alpha_composite(molecules, circle)
    image = image.crop(image.getbbox())
    image.quantize(colors=256, method=Image.Quantize.FASTOCTREE, dither=Image.Dither.NONE).save(OUT, optimize=True)
    print(f"wrote {OUT} ({image.width} x {image.height} px, {OUT.stat().st_size / 1e6:.2f} MB)")


if __name__ == "__main__":
    main()
