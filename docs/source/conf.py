from __future__ import annotations

import os
import subprocess
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]

project = "Q and QligFEP"
author = "The Q and QligFEP developers"
copyright = "The Q and QligFEP developers"


def repository_version() -> str:
    """Return a useful documentation version without importing QligFEP."""
    if os.environ.get("READTHEDOCS_VERSION_NAME"):
        return os.environ["READTHEDOCS_VERSION_NAME"]
    try:
        return subprocess.check_output(
            ["git", "describe", "--tags", "--always", "--dirty"],
            cwd=ROOT,
            text=True,
            stderr=subprocess.DEVNULL,
        ).strip()
    except (OSError, subprocess.CalledProcessError):
        return "development"


release = repository_version()
version = release

extensions = [
    "myst_parser",
    "sphinx_copybutton",
    "sphinx_design",
    "sphinxcontrib.mermaid",
    "sphinx.ext.extlinks",
]

source_suffix = {".md": "markdown"}
master_doc = "index"
exclude_patterns = ["_build", "Thumbs.db", ".DS_Store"]

myst_enable_extensions = [
    "colon_fence",
    "deflist",
    "dollarmath",
    "fieldlist",
    "substitution",
]
myst_heading_anchors = 4
myst_url_schemes = ("http", "https", "mailto")

html_theme = "furo"
html_title = "Q and QligFEP"
html_static_path = ["_static"]
html_css_files = ["custom.css"]
html_theme_options = {
    "source_repository": "https://github.com/qusers/Q/",
    "source_branch": "main",
    "source_directory": "docs/source/",
}

extlinks = {
    "repo": ("https://github.com/qusers/Q/blob/main/%s", "%s"),
    "issue": ("https://github.com/qusers/Q/issues/%s", "issue #%s"),
}

copybutton_prompt_text = r">>> |\.\.\. |\$ |In \[\d*\]: | {2,5}\.\.\.: | {5,8}: "
copybutton_prompt_is_regexp = True

# ACS blocks automated HEAD/GET requests with HTTP 403 even though DOI links are
# valid in a browser. Keep these citations but exclude that bot-specific failure.
linkcheck_ignore = [r"https://doi\.org/10\.1021/.*"]
