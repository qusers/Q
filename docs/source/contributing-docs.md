# Documentation and Read the Docs

## Local build

```bash
python -m pip install -r docs/requirements.txt
sphinx-build -W --keep-going -b html docs/source docs/build/html
```

The warning-as-error build is the required local check. Run the external link
checker separately:

```bash
sphinx-build -b linkcheck docs/source docs/build/linkcheck
```

## Source layout

- `docs/source/` contains the Read the Docs navigation and topic pages.
- `docs/source/concepts/` explains how Q does FEP.
- `docs/source/manual/` is the Q manual, one page per chapter.
- `docs/scripts/` holds the scripts that generate figures.
- `.readthedocs.yml` defines the hosted build.

Keep generated HTML under `docs/build/` out of version control.

## Publish on Read the Docs

1. Sign in to [Read the Docs](https://readthedocs.org/).
2. Import `qusers/Q` (or the final repository URL).
3. Keep the configuration-file path as `.readthedocs.yml`.
4. Build the default branch.
5. In the build log, confirm Python 3.13, `docs/requirements.txt`, and
   `docs/source/conf.py` were used.
6. Configure the desired public project slug and default version in the Read the
   Docs project settings.

The repository configuration is sufficient for a documentation build; it does
not install QligFEP or compile the simulation engine on the documentation worker.
This keeps hosted builds independent of heavy scientific runtime dependencies.

## Authoring rules

- Document behavior supported by the current source, not historical plans.
- Put exact option/default changes in the CLI reference in the same change.
- Include units and define edge direction for scientific values.
- Distinguish setup, submission, simulation, and analysis steps.
- Label destructive commands and experimental methods.
- Use relative links for repository documentation.
- Do not commit generated build output.
