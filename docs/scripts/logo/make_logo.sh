#!/usr/bin/env bash
# Rebuild the chosen Q logo (cdk2 ligand 1h1q, polar hydrogens only, in its water sphere).
# Run from the repository root:  bash docs/scripts/logo/make_logo.sh
# Writes docs/source/_static/q-logo.svg, the logo shown above the sidebar of the documentation.
# Requires: pymol (with numpy), and for --retrace also ImageMagick + potrace.
set -euo pipefail
HERE=docs/scripts/logo
OUT=docs/source/_static/q-logo.svg

if [[ "${1:-}" == "--retrace" ]]; then
    # 1) Trace the red Q and its bowl from the original raster logo (3x upsampled).
    tmp=$(mktemp -d)
    magick "$HERE/original-q-logo.png" -alpha off \
        -fx '(r>0.6 && g<0.35 && b<0.35) ? 0 : 1' -morphology Close Disk:2 \
        -resize 300% -blur 0x2 -threshold 50% "$tmp/red.pbm"
    # drop the red water oxygens of the original picture (small connected components)
    magick "$tmp/red.pbm" -define connected-components:area-threshold=20000 \
        -define connected-components:mean-color=true -connected-components 8 \
        -threshold 50% "$tmp/red-clean.pbm"
    potrace "$tmp/red-clean.pbm" -s -t 2000 -a 1.2 -O 0.4 --flat -o "$HERE/q-outline-traced.svg"
    # bowl = interior of the Q (background flood-filled away)
    magick "$tmp/red-clean.pbm" -negate -fill white -draw "color 0,0 floodfill" "$tmp/bowl.pbm"
    potrace "$tmp/bowl.pbm" -s -t 2000 -a 1.2 -O 0.4 --flat -o "$HERE/q-bowl-traced.svg"
fi

# 2) Ligand in the 20 Å water sphere of test/q6, projected to SVG.
pymol -cq "$HERE/build_sphere_logo.py" -- \
    test/qligfep/resources/cdk2_ligands.sdf 1h1q \
    "$OUT" \
    scale=74 sphere=20 stick=0.17 keep=1 slab=5 halo=2.0 wrot=37 rot=-20 hyd=polar
