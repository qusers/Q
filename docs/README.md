# Q documentation

The documentation is built with Sphinx, MyST Markdown, and the Furo theme. It
covers the Q engine and the QligFEP and QresFEP workflows.

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
