"""Figure for docs/source/concepts/spherical-boundary.md: periodic system and Q sphere.

Places two existing PyMOL renders of the adenosine A2A receptor side by side:
the full periodic membrane system and the sphere that Q simulates. The renders
are read in place from a folder that is not part of the repository; nothing is
re-rendered.

    python docs/scripts/compose_sphere_vs_box.py --memprep-dir /path/to/GPCR-A2A-Deflorian/memprep
"""

from __future__ import annotations

import argparse
from pathlib import Path

from PIL import Image

OUT = Path(__file__).resolve().parents[1] / "source" / "_static" / "figures" / "sphere-vs-box.png"
PANELS = ("membrane-protein.png", "simulated-prot-part.png")
HEIGHT = 640  # px
GAP = 60  # px


def load(path: Path) -> Image.Image:
    """Crop to the rendered content and scale to the panel height (premultiplied alpha)."""
    image = Image.open(path).convert("RGBA")
    image = image.crop(image.getbbox()).convert("RGBa")
    width = round(image.width * HEIGHT / image.height)
    return image.resize((width, HEIGHT), Image.LANCZOS).convert("RGBA")


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--memprep-dir", type=Path, required=True)
    args = parser.parse_args()
    panels = [load(args.memprep_dir / name) for name in PANELS]
    canvas = Image.new("RGBA", (sum(p.width for p in panels) + GAP * (len(panels) - 1), HEIGHT), (255, 255, 255, 0))
    x = 0
    for panel in panels:
        canvas.paste(panel, (x, 0))
        x += panel.width + GAP
    OUT.parent.mkdir(parents=True, exist_ok=True)
    # 256 colours with alpha keep the file below 1 MB; the renders are not photographs
    canvas.quantize(colors=256, method=Image.Quantize.FASTOCTREE, dither=Image.Dither.NONE).save(OUT, optimize=True)
    print(f"wrote {OUT} ({canvas.width} x {canvas.height} px, {OUT.stat().st_size / 1e6:.2f} MB)")


if __name__ == "__main__":
    main()
