# Q and QligFEP documentation

The documentation is built with Sphinx, MyST Markdown, and the Furo theme. The
maintained engine manual remains at
`docs/source/q-engine/manual.md` and is included in this site.

## Build locally

From the repository root:

```bash
python -m pip install -r docs/requirements.txt
sphinx-build -W --keep-going -b html docs/source docs/build/html
```

Or from `docs/`:

```bash
make html SPHINXOPTS="-W --keep-going"
```

Open `docs/build/html/index.html`. Run `make linkcheck` separately to test
external links. Generated files under `docs/build/` should not be committed.
