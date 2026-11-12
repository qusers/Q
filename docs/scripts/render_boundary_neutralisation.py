"""Figure for docs/source/concepts/electrostatics.md: ionisable residues by treatment.

Renders the Tyk2 structure of the benchmark set with PyMOL, as prepared by
``qprep_prot``: a whole-protein cartoon for context, and the side chains of
the ionisable residues whose ionisable group lies within 6 A of the plane
through the sphere centre. The residues that were neutralised are read from
``preparation.json``; which atoms are excluded comes from the qprep mask (see
``render_sphere_regions.py`` for the qprep input that writes it). The
structures are read in place and are not part of the repository.

The image is built from separately rendered layers in one orthoscopic view
(cartoon, surface of the excluded atoms in the slab, side chains, ligand,
circle), so the context cannot hide the residues. The dashed circle lies in
the plane through the sphere centre perpendicular to the view.

    python docs/scripts/render_boundary_neutralisation.py \
        --benchmark-dir /path/to/perturbations/tyk2 --free /path/to/notexcluded.pdb
"""
from __future__ import annotations

import argparse
import json
import math
import tempfile
from pathlib import Path

import numpy as np
from PIL import Image, ImageChops, ImageFilter
from pymol import cgo, cmd

OUT = Path(__file__).resolve().parents[1] / "source" / "_static" / "figures" / "boundary-neutralisation.png"
WORK = Path(tempfile.gettempdir()) / "boundary-neutralisation-layers"  # intermediate renders
SIZE = 1500
NORMAL = (1., 0., 0.)
UP = (0., 0., 1.)
COLORS = {"acid": (0.86, 0.25, 0.29), "base": (0.18, 0.55, 0.82),
          "neutral": (0.88, 0.83, 0.71), "backbone": (0.65, 0.70, 0.73),
          "outside": (0.80, 0.77, 0.69), "inside": (0.60, 0.72, 0.77),
          "ligand_color": (0.98, 0.49, 0.07)}
# Heavy atoms of the terminal ionisable moiety (also used for neutral forms).
GROUPS = {"ASP": "CG OD1 OD2", "GLU": "CD OE1 OE2", "LYS": "NZ",
          "ARG": "NE CZ NH1 NH2", "HIP": "CG ND1 CD2 CE1 NE2"}


def select_indices(name, indices):
    cmd.select_list(name, "protein", indices, mode="index")


def circle(radius):
    obj = []
    color = (0.46, 0.48, 0.51)
    for i in range(64):
        angles = np.linspace(2 * math.pi * i / 64, 2 * math.pi * (i + .57) / 64, 5)
        for a, b in zip(angles[:-1], angles[1:]):
            p = (radius * math.cos(a), radius * math.sin(a), 0.)
            q = (radius * math.cos(b), radius * math.sin(b), 0.)
            obj += [cgo.CYLINDER, *p, *q, .080, *color, *color]
    return obj


def render(tag, opacity=1., outline=.75):
    """Keep intermediate rays local; add only an exterior 1 px silhouette."""
    path = WORK / (tag + ".png")
    cmd.png(str(path), width=SIZE, height=SIZE, ray=1, dpi=200)
    im = Image.open(path).convert("RGBA")
    alpha = im.getchannel("A")
    ring = ImageChops.subtract(alpha.filter(ImageFilter.MaxFilter(3)), alpha)
    edge = Image.new("RGBA", im.size, (62, 67, 72, 0))
    edge.putalpha(ring.point(lambda x: round(x * outline)))
    im = Image.alpha_composite(edge, im)
    im.putalpha(im.getchannel("A").point(lambda x: round(x * opacity)))
    return im


def save(image):
    box = image.getbbox()
    pad = 22
    box = (max(0, box[0]-pad), max(0, box[1]-pad),
           min(SIZE, box[2]+pad), min(SIZE, box[3]+pad))
    image = image.crop(box)
    image = image.resize(tuple(round(v * 1450 / max(image.size)) for v in image.size), Image.Resampling.LANCZOS)
    # Indexed RGBA keeps this flat-shaded illustration compact, including alpha.
    image.quantize(colors=256, method=Image.Quantize.FASTOCTREE,
                   dither=Image.Dither.NONE).save(OUT, optimize=True)
    print(f"Wrote {OUT.name}: {image.size}, {OUT.stat().st_size / 1000:.0f} kB", flush=True)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--benchmark-dir", type=Path, required=True)
    parser.add_argument("--free", type=Path, required=True)
    args = parser.parse_args()
    WORK.mkdir(parents=True, exist_ok=True)
    prep = json.loads((args.benchmark_dir / "preparation.json").read_text())
    centre = np.asarray(prep["legacy_cog"])
    radius = float(prep["nominal_radius"])
    changes = prep["neutralization"]["modifications"]
    cmd.reinitialize()
    cmd.set("auto_zoom", 0)
    cmd.load(str(args.benchmark_dir / "protein.pdb"), "protein")
    cmd.load(str(args.benchmark_dir / "ejm_31.pdb"), "ligand")
    cmd.load(str(args.free), "mask")
    cmd.remove("hydrogens")
    # Coordinate matching handles qprep's absent chain IDs and reordered atoms.
    cmd.select("free_atoms", "protein and (mask around 0.01)")
    cmd.select("excluded_atoms", "protein and not free_atoms")
    normal = np.asarray(NORMAL, dtype=float)
    normal /= np.linalg.norm(normal)
    right = np.cross(UP, normal)
    right /= np.linalg.norm(right)
    rotation = np.array([right, np.cross(normal, right), normal])
    for obj in ("protein", "ligand", "mask"):
        cmd.load_coords((cmd.get_coords(obj) - centre) @ rotation.T, obj)
    residues = {}
    for atom in cmd.get_model("protein").atom:
        residues.setdefault(f"{atom.chain}:{atom.resi}", []).append(atom)
    assert set(changes) <= set(residues), "Preparation residue missing from protein"
    classes = {"acid": [], "base": [], "neutral": []}
    moieties, shown = [], {k: [] for k in classes}
    for key, atoms in residues.items():
        resn = atoms[0].resn
        if key in changes:
            assert resn in (changes[key]["original"], changes[key]["modified"]), key
            category, original = "neutral", changes[key]["original"]
        elif resn in GROUPS:
            category, original = ("acid" if resn in ("ASP", "GLU") else "base"), resn
        else:
            continue
        group = [a for a in atoms if a.name in GROUPS[original].split()]
        assert len(group) == len(GROUPS[original].split()), (key, original)
        classes[category].append(key)
        # Visibility selection only; this never determines the charge class.
        if abs(np.mean([a.coord[2] for a in group])) > 6.:
            continue
        shown[category].append(key)
        moieties.extend(a.index for a in group)
    assert [len(classes[k]) for k in ("acid", "base", "neutral")] == [23, 20, 32]
    for category, keys in shown.items():
        select_indices(category, [a.index for key in keys for a in residues[key]])
    select_indices("moieties", moieties)
    print("Residues drawn: " + str({k: len(v) for k, v in shown.items()}), flush=True)
    for name, rgb in COLORS.items():
        cmd.set_color(name, rgb)
    cmd.hide("everything")
    cmd.dss("protein")
    cmd.color("backbone", "protein")
    cmd.color("outside", "excluded_atoms")
    cmd.color("ligand_color", "ligand")
    cmd.set("cartoon_fancy_helices", 1)
    cmd.set("cartoon_flat_sheets", 1)
    cmd.set("cartoon_loop_radius", .18)
    cmd.set("cartoon_oval_width", .22)
    cmd.set("cartoon_oval_length", 1.05)
    cmd.set("cartoon_rect_width", .22)
    cmd.set("cartoon_rect_length", 1.1)
    cmd.set("stick_radius", .19)
    cmd.set("stick_ball", 0)
    cmd.set("valence", 0)
    cmd.set("orthoscopic", 1)
    cmd.set("depth_cue", 0)
    cmd.set("ray_shadows", 0)
    cmd.set("specular", .05)
    cmd.set("ambient", .65)
    cmd.set("direct", .35)
    cmd.set("ray_opaque_background", 0)
    cmd.set("antialias", 2)
    cmd.set("surface_quality", 1)
    cmd.set("max_threads", 4)
    xy = cmd.get_coords("protein or ligand")[:, :2]
    low = np.minimum(xy.min(0)-3, -radius-2)
    high = np.maximum(xy.max(0)+3, radius+2)
    mid = (low+high)/2
    cmd.set_view([1,0,0, 0,1,0, 0,0,1, 0,0,-200, *mid,0, 100,300,0])
    cmd.pseudoatom("frame", pos=[*mid, 0.])
    cmd.zoom("frame", float(max(high-low)/2))
    cmd.hide("everything", "frame")
    view = cmd.get_view()
    # Separate layers keep the protein context from obscuring selected residues.
    cmd.show("cartoon", "protein")
    context = render("context", opacity=.34, outline=.25)
    cmd.hide("everything")
    cmd.create("outer_slice", "excluded_atoms and z > -6 and z < 6")
    cmd.show("surface", "outer_slice")
    context = Image.alpha_composite(context, render("slab_surface", opacity=.20, outline=.3))
    cmd.hide("everything")
    for category in classes:
        cmd.color(category, category)
    cmd.show("sticks", "(acid or base or neutral) and (sidechain or name CA)")
    highlights = render("highlights")
    result = Image.alpha_composite(context, highlights)
    cmd.hide("everything")
    cmd.show("sticks", "ligand")
    cmd.set("stick_radius", .13, "ligand")
    result = Image.alpha_composite(result, render("ligand"))
    cmd.hide("everything")
    cmd.load_cgo(circle(radius), "boundary_circle")
    cmd.set_view(view)
    result = Image.alpha_composite(result, render("circle", outline=0))
    save(result)


if __name__ == "__main__":
    main()
