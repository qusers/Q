"""Build a fully-vector Q logo: a ligand in a spherical water boundary inside the Q.

PyMOL is used to assemble and orient the system (ligand dropped in the centre of the
pre-equilibrated water sphere from test/q6, clashing waters removed). Atom positions are
then projected to camera space and drawn as SVG primitives (shaded circles for the
ligand, two-tone sticks for the waters), painted back-to-front.

Run from the repository root:
    pymol -cq docs/scripts/logo/build_sphere_logo.py -- <sdf> <title> <out.svg> [key=value ...]

The recipe of the logo used in the documentation is in make_logo.sh.
"""

import re
import sys
from pathlib import Path

import numpy as np
from pymol import cmd

ROOT = Path.cwd()
ASSETS = ROOT / "docs/scripts/logo"
WATERS = ROOT / "test/q6/test1/prep/lig_w.pdb"

BOWL_BBOX = (339, 294, 2190, 2277)  # bounding box of the traced bowl (logo px)

ELEM = {  # fill, highlight, rim, vdW radius (Å)
    "C": ("#4a4e55", "#8a9099", "#26292d", 1.70),
    "N": ("#2f56d6", "#7d9cff", "#1a2f80", 1.55),
    "O": ("#d92a2a", "#ff8080", "#7a1010", 1.52),
    "S": ("#d6b21f", "#ffe680", "#7a6410", 1.80),
    "Cl": ("#2fa84f", "#8be0a2", "#17592a", 1.75),
    "F": ("#5fbf8f", "#b5f0d2", "#2e6b4c", 1.47),
    "H": ("#e9ecef", "#ffffff", "#9aa1a8", 1.20),
}


def parse_args(argv):
    sdf, title, out, *rest = argv
    opts = dict(
        sphere=20.0,  # sphere radius (Å); 20 Å is the Q water sphere of test/q6
        slab=99.0,  # keep waters with |z| < slab (Å)
        clash=2.6,  # remove waters within this distance of the ligand (Å)
        scale=62.0,  # px per Å
        rot=0.0,  # in-plane rotation (deg)
        tilt=0.0,  # rotation about x (deg)
        spin=0.0,  # rotation about y (deg)
        cx=1290.0,  # sphere centre in logo px
        cy=1255.0,
        stick=0.30,  # water stick radius (Å)
        shell=0,  # draw an explicit sphere disc/shell (off: the Q bowl is the boundary)
        keep=1.0,  # fraction of waters kept (random thinning)
        edgefade=0.0,  # extra fade towards the sphere edge (0: uniform)
        halo=3.2,
        hole=0.0,  # also drop waters behind the ligand within this xy distance (Å)  # no waters in front of the ligand within this xy distance (Å)
        seed=3,
        hyd="all",  # ligand hydrogens: all | polar (hide H on carbon) | none
        style="spheres",  # ligand representation: spheres | sticks
        lstick=0.30,  # ligand stick radius (Å)
        ltilt=0.0,  # rotate the ligand only about the camera x axis (deg)
        lzoom=1.0,
        ghost="",  # title of a second ligand (same SDF) overlaid as a transparent ghost
        c1="#5b6470",  # carbon colour of the main ligand in ghost mode
        c2="#e0a43a",  # carbon colour of the ghost ligand
        gop=0.78,  # opacity of each ligand in ghost mode
        shine=0,
        gx=0.0,  # screen-space offset of the ghost ligand (Å, before zoom)
        gy=0.0,  # highlight stripe along the sticks  # magnify the ligand only (keeps the water look unchanged)
        wrot=37.0,  # rotation of the water sphere (deg) about waxis
        waxis=[0.31, 0.83, 0.47],
    )
    for kv in rest:
        k, v = kv.split("=")
        opts[k] = type(opts[k])(float(v)) if isinstance(opts[k], (int, float)) else v
    assert opts["hyd"] in ("all", "polar", "none") and opts["style"] in ("spheres", "sticks")
    return Path(sdf), title, Path(out), opts


def build_system(sdf, title, o):
    cmd.reinitialize()
    cmd.load(str(WATERS), "wsrc")
    cmd.alter("wsrc", "elem = name[0]")
    centre = cmd.get_coords("wsrc and resn HOH and name O*").mean(0)
    cmd.create("wat", "wsrc and resn HOH")
    cmd.create("oldlig", "wsrc and resn LIG")
    cmd.delete("wsrc")
    fill_old_cavity(centre)
    # the sphere is cut from a replicated box: rotate it so no lattice axis faces the
    # camera (otherwise the waters line up in visible rows)
    cmd.rotate(o["waxis"], o["wrot"], "wat", camera=0, origin=list(map(float, centre)))

    cmd.load(str(sdf), "lib", discrete=1)
    for state in range(1, cmd.count_states("lib") + 1):
        if cmd.get_title("lib", state) == title:
            cmd.create("lig", "lib", state, 1)
            break
    else:
        raise SystemExit(f"title {title!r} not found in {sdf}")
    if o["ghost"]:
        for state in range(1, cmd.count_states("lib") + 1):
            if cmd.get_title("lib", state) == o["ghost"]:
                cmd.create("lig2", "lib", state, 1)
                break
        else:
            raise SystemExit(f"ghost {o['ghost']!r} not found in {sdf}")
    cmd.delete("lib")
    # both ligands share the docking frame, so move them together
    shift = centre - cmd.get_coords("lig").mean(0)
    cmd.translate(list(map(float, shift)), "lig*", camera=0)
    cmd.remove(f"byres (wat within {o['clash']} of (lig* and not elem H))")

    cmd.orient("lig")
    cmd.turn("x", o["tilt"])
    cmd.turn("y", o["spin"])
    cmd.turn("z", o["rot"])
    # tilt only the ligand (camera axes, about its centroid): reveals rings without
    # changing how the water lattice is viewed
    if o["ltilt"]:
        cen = list(map(float, cmd.get_coords("lig").mean(0)))
        cmd.rotate("x", o["ltilt"], "lig*", camera=1, origin=cen)
    # hide hydrogens only after orienting, so the view doesn't depend on them
    if o["hyd"] == "polar":
        cmd.remove("lig* and elem H and (neighbor elem C)")
    elif o["hyd"] == "none":
        cmd.remove("lig* and elem H")
    return centre


def fill_old_cavity(centre, min_oo=2.7):
    """Fill the cavity left by the ligand the sphere was built around.

    A copy of the sphere, rotated by an arbitrary angle about the centre (180° maps the
    periodic water lattice onto itself), donates the waters whose O falls inside the
    cavity and doesn't clash with existing waters.
    """
    cmd.create("donor", "wat")
    cmd.alter("donor", "resi = str(int(resi) + 100000)")
    c = list(map(float, centre))
    cmd.rotate([1, 2, 3], 73, "donor", camera=0, origin=c)
    cmd.select("incav", "byres (donor and name O* and (donor within 4.0 of oldlig))")
    cmd.remove("donor and not incav")
    cmd.delete("incav")
    cmd.remove(f"byres (donor and name O* and (donor within {min_oo} of (wat and name O*)))")
    # greedy de-clash among the donated waters themselves
    kept = []
    for resi, xyz in sorted(
        (a.resi, a.coord) for a in cmd.get_model("donor and name O*").atom
    ):
        if all(np.linalg.norm(np.subtract(xyz, k)) >= min_oo for k in kept):
            kept.append(xyz)
        else:
            cmd.remove(f"donor and resi {resi}")
    print(f"cavity filled with {len(kept)} waters")
    cmd.create("wat", "wat or donor")
    cmd.delete("donor")
    cmd.delete("oldlig")


def camera_coords(sel):
    view = cmd.get_view()
    rot = np.array(view[:9]).reshape(3, 3)
    origin = np.array(view[12:15])
    xyz = cmd.get_coords(sel)
    return (xyz - origin) @ rot  # x right, y up, z towards viewer


def atoms(sel):
    model = cmd.get_model(sel)
    return [(a.symbol or a.name[0], a.resi, a.name, a.index) for a in model.atom]


def mix(hex_a, hex_b, t):
    a = np.array([int(hex_a[i : i + 2], 16) for i in (1, 3, 5)])
    b = np.array([int(hex_b[i : i + 2], 16) for i in (1, 3, 5)])
    c = np.round(a + (b - a) * t).astype(int)
    return "#" + "".join(f"{v:02x}" for v in c)


def path_d(svg_file):
    text = Path(svg_file).read_text()
    d = re.search(r'<path d="([^"]+)"', text, re.S).group(1)
    return " ".join(d.split())


def ligand_sticks(sel, lig_at, lig_xyz, px, w, carbon=None, shine=False):
    """Licorice sticks: two-tone half-bonds, each with a rim in a darker shade of itself."""
    bonds = [b.index for b in cmd.get_model(sel).bond]
    assert len(lig_at) == len(lig_xyz)
    elems = [el if el in ELEM else "C" for el, *_ in lig_at]

    def colour(el):
        return carbon if (carbon and el == "C") else ELEM[el][0]

    segs = []  # (z, (x0, y0), (x1, y1), element)
    for i, j in bonds:
        mid = (lig_xyz[i] + lig_xyz[j]) / 2
        for a in (i, j):
            segs.append(((lig_xyz[a][2] + mid[2]) / 2, px(lig_xyz[a]), px(mid), elems[a]))
    segs.sort(key=lambda t: t[0])
    line = lambda a, b, dx=0.0: f"M{a[0] + dx:.1f} {a[1] + dx:.1f}L{b[0] + dx:.1f} {b[1] + dx:.1f}"
    rim = "".join(f'<path d="{line(a, b)}" stroke="{mix(colour(el), "#000000", 0.4)}"/>' for _, a, b, el in segs)
    body = "".join(f'<path d="{line(a, b)}" stroke="{colour(el)}"/>' for _, a, b, el in segs)
    out = f'<g stroke-width="{w + 12:.1f}">{rim}</g><g stroke-width="{w:.1f}">{body}</g>'
    if shine:
        hi = "".join(
            f'<path d="{line(a, b, -0.18 * w)}" stroke="{ELEM[el][1]}"/>' for _, a, b, el in segs
        )
        out += f'<g stroke-width="{0.3 * w:.1f}" stroke-opacity="0.55">{hi}</g>'
    return out


def build_svg(centre, o):
    # Camera-space coordinates, centred on the sphere centre.
    cam_c = (centre - np.array(cmd.get_view()[12:15])) @ np.array(cmd.get_view()[:9]).reshape(3, 3)
    lig_xyz = camera_coords("lig") - cam_c
    lig_at = atoms("lig")
    # optional artistic magnification of the ligand about its centroid (waters unchanged)
    zc = lig_xyz.mean(0)
    lig_xyz = zc + (lig_xyz - zc) * o["lzoom"]
    ghost = bool(o["ghost"])
    if ghost:
        lig2_at = atoms("lig2")
        lig2_xyz = zc + (camera_coords("lig2") - cam_c - zc) * o["lzoom"]
        # benchmark FEP poses share an identical core: an optional small screen-space
        # offset makes the two "ghosts" read as separate molecules
        lig2_xyz = lig2_xyz + np.array([o["gx"], o["gy"], -0.5])
    wat_xyz = camera_coords("wat") - cam_c
    wat_at = atoms("wat")

    rng = np.random.default_rng(int(o["seed"]))
    lig_heavy = np.array([p for (el, *_), p in zip(lig_at, lig_xyz) if el != "H"])
    if ghost:
        lig_heavy = np.vstack([lig_heavy, [p for (el, *_), p in zip(lig2_at, lig2_xyz) if el != "H"]])
    # group waters by residue
    waters = {}
    for (el, resi, name, _), p in zip(wat_at, wat_xyz):
        waters.setdefault(resi, []).append((el, p))
    R, S = o["sphere"], o["scale"]

    def px(p):
        return o["cx"] + p[0] * S, o["cy"] - p[1] * S

    items = []  # (depth, svg)
    bg = "#ffffff"
    for resi, mol in waters.items():
        ox = [p for el, p in mol if el == "O"][0]
        r = np.linalg.norm(ox)
        if r > R or abs(ox[2]) > o["slab"] or rng.random() > o["keep"]:
            continue
        # keep the ligand unobstructed: drop waters in front of it in projection
        d2 = np.linalg.norm(lig_heavy[:, :2] - ox[:2], axis=1)
        if np.any((d2 < o["halo"]) & (ox[2] > lig_heavy[:, 2] - 1.0)):
            continue
        # ...and directly behind it (visible through transparent ghosts)
        if np.any(d2 < o["hole"]):
            continue
        # fade: back waters and waters near the shell get lighter
        depth_t = np.clip((R - ox[2]) / (2 * R), 0, 1)  # 0 front, 1 back
        edge_t = np.clip((r - 0.7 * R) / (0.3 * R), 0, 1)
        fade = 0.05 + 0.45 * depth_t + o["edgefade"] * edge_t
        fade = min(fade, 0.9)
        x0, y0 = px(ox)
        if not (BOWL_BBOX[0] - 60 < x0 < BOWL_BBOX[2] + 60 and BOWL_BBOX[1] - 60 < y0 < BOWL_BBOX[3] + 60):
            continue  # would be clipped away by the bowl anyway
        w = o["stick"] * 2 * S
        parts = []
        for el, h in mol:
            if el != "H":
                continue
            mid = (ox + h) / 2
            xm, ym = px(mid)
            xh, yh = px(h)
            edge = mix("#5a5f66", bg, fade)
            parts.append(
                f'<path d="M{x0:.1f} {y0:.1f}L{xh:.1f} {yh:.1f}" stroke="{edge}" stroke-width="{w + 6:.1f}"/>'
            )
        for el, h in mol:
            if el != "H":
                continue
            mid = (ox + h) / 2
            xm, ym = px(mid)
            xh, yh = px(h)
            parts.append(
                f'<path d="M{xm:.1f} {ym:.1f}L{xh:.1f} {yh:.1f}" stroke="{mix("#f2f3f5", bg, fade)}" stroke-width="{w:.1f}"/>'
                f'<path d="M{x0:.1f} {y0:.1f}L{xm:.1f} {ym:.1f}" stroke="{mix("#d8262b", bg, fade)}" stroke-width="{w:.1f}"/>'
            )
        items.append((ox[2], "".join(parts)))

    if o["style"] == "spheres":
        for (el, *_), p in zip(lig_at, lig_xyz):
            el = el if el in ELEM else "C"
            fill, hi, rim, rad = ELEM[el]
            x, y = px(p)
            rr = rad * S
            items.append(
                (p[2], f'<circle cx="{x:.1f}" cy="{y:.1f}" r="{rr:.1f}" fill="url(#g{el})" stroke="{rim}" stroke-width="5"/>')
            )
    else:
        w = o["lstick"] * 2 * S
        if ghost:
            zmax = max(lig_xyz[:, 2].max(), lig2_xyz[:, 2].max())
            g1 = ligand_sticks("lig", lig_at, lig_xyz, px, w, o["c1"], o["shine"])
            g2 = ligand_sticks("lig2", lig2_at, lig2_xyz, px, w, o["c2"], o["shine"])
            items.append((zmax, f'<g opacity="{o["gop"]}">{g1}</g><g opacity="{o["gop"]}">{g2}</g>'))
        else:
            items.append((lig_xyz[:, 2].max(), ligand_sticks("lig", lig_at, lig_xyz, px, w, shine=o["shine"])))

    items.sort(key=lambda t: t[0])
    grads = "".join(
        f'<radialGradient id="g{el}" cx="0.36" cy="0.32" r="0.75">'
        f'<stop offset="0" stop-color="{hi}"/><stop offset="0.45" stop-color="{fill}"/>'
        f'<stop offset="1" stop-color="{rim}"/></radialGradient>'
        for el, (fill, hi, rim, _) in ELEM.items()
    )
    q_d = path_d(ASSETS / "q-outline-traced.svg")
    bowl_d = path_d(ASSETS / "q-bowl-traced.svg")
    potrace_tf = "translate(0,3096) scale(0.1,-0.1)"
    sphere_r = R * S
    shell = ""
    if o["shell"]:
        shell = (
            f'<circle cx="{o["cx"]}" cy="{o["cy"]}" r="{sphere_r:.1f}" fill="url(#sph)"/>'
            f'<circle cx="{o["cx"]}" cy="{o["cy"]}" r="{sphere_r:.1f}" fill="none" stroke="#9fb3c8" stroke-width="8" stroke-dasharray="2 26" stroke-linecap="round"/>'
        )
    return f"""<svg xmlns="http://www.w3.org/2000/svg" viewBox="0 0 2550 3096">
<defs>{grads}
<radialGradient id="sph" cx="0.4" cy="0.38" r="0.7"><stop offset="0" stop-color="#ffffff"/><stop offset="0.75" stop-color="#eef3f8"/><stop offset="1" stop-color="#dbe5ef"/></radialGradient>
<clipPath id="bowl"><path transform="{potrace_tf}" d="{bowl_d}"/></clipPath>
</defs>
<path transform="{potrace_tf}" d="{bowl_d}" fill="#ffffff"/>
<g clip-path="url(#bowl)" stroke-linecap="round" fill="none">
{shell}
{"".join(s for _, s in items)}
</g>
<path transform="{potrace_tf}" d="{q_d}" fill="#e3061b" stroke="#5a0008" stroke-width="40"/>
</svg>
"""


def main():
    sdf, title, out, o = parse_args(sys.argv[1:])
    centre = build_system(sdf, title, o)
    out.write_text(build_svg(centre, o))
    print(f"wrote {out}")


main()
