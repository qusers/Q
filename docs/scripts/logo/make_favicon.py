"""Write the favicon of the documentation: the red Q of the logo, without the sphere inside.

Uses the traced outline next to this script and needs only the standard library.
Run from the repository root (make_logo.sh does this and also writes the PNG):
    python docs/scripts/logo/make_favicon.py
"""

import re
from pathlib import Path

HERE = Path(__file__).resolve().parent
OUT = HERE.parents[1] / "source" / "_static" / "q-favicon.svg"
WIDTH, HEIGHT = 2550, 3096  # size of the traced logo
BOLD = 90  # px of the traced logo added to the stroke of the Q, so it stays visible at 16 px
SIDE = HEIGHT + BOLD  # the favicon is a centred square

path = " ".join(re.search(r'<path d="([^"]+)"', (HERE / "q-outline-traced.svg").read_text(), re.S).group(1).split())
OUT.write_text(
    f'<svg xmlns="http://www.w3.org/2000/svg" viewBox="{-(SIDE - WIDTH) // 2} {-BOLD // 2} {SIDE} {SIDE}">\n'
    f'<path transform="translate(0,{HEIGHT}) scale(0.1,-0.1)" d="{path}" fill="#e3061b" stroke="#e3061b" '
    f'stroke-width="{BOLD * 10}" stroke-linejoin="round"/>\n</svg>\n'
)
print(f"wrote {OUT}")
