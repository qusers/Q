"""Light and dark variants of a figure from one ``figure`` or ``image`` directive.

When a page includes ``<name>-light.svg`` and ``<name>-dark.svg`` exists next
to it, the dark file is added after the light one. Furo's ``only-light`` and
``only-dark`` classes then show the file that matches the active theme, so the
caption is written once.
"""

from __future__ import annotations

from pathlib import Path

from docutils import nodes
from sphinx.application import Sphinx
from sphinx.transforms import SphinxTransform

LIGHT, DARK = "-light.", "-dark."


class ThemedFigures(SphinxTransform):
    default_priority = 400

    def apply(self, **kwargs) -> None:
        source_dir = Path(self.env.doc2path(self.env.docname)).parent
        for image in list(self.document.findall(nodes.image)):
            uri = image["uri"]
            if LIGHT not in uri or "only-light" in image["classes"]:
                continue
            dark_uri = DARK.join(uri.rsplit(LIGHT, 1))
            if not (source_dir / dark_uri).is_file():
                continue
            dark = image.deepcopy()
            dark["uri"] = dark_uri
            dark["classes"].append("only-dark")
            image["classes"].append("only-light")
            image.parent.insert(image.parent.index(image) + 1, dark)


def setup(app: Sphinx) -> dict:
    app.add_transform(ThemedFigures)
    return {"parallel_read_safe": True, "parallel_write_safe": True}
