# CHANGELOG - QligFEP refactoring:
## New packaging structure:
```text
Q/
├── LICENSE
├── pyproject.toml
├── setup.py # needed for compiling Q
├── README.md
├── src/
│   ├── Q6/
│   │   ├── <filename>.f90
│   │   └── makefile
│   ├── QligFEP/ # old name: qligfep-old /
│   │   ├── __init__.py
│   │   ├── <other modules/subpackages>
│   │   └── main.py # command line interface # used to be under bin/
│   ├── Qgpu/ # old name: share/
│   │   ├── __init__.py
│   │   └── other_stuff.py # other modules within package
│   ├── core/
│   │   ├── <modules>.cu
│   │   └── <modules>.h
│   └── qligfep-newbin-unfinished/ **
│       ├── __init__.py
│       └── other_stuff.py # other modules within package
└── tests/
```
** unfinished re-implementation of qligfep by willem. Used to be under `src/`.

## Changes
- Moved around files to follow the [new packaging structure](#new-packaging-structure).
- Updated `.gitignore` with templates for `python` and for `fortran` to avoid adding unnecessary files to Git.
- Removed all files/directories such as `.egg-info` and `.ipynb_checkpoints`.
- Changed all `sys.path.insert` statements to the correct respective import statement:
  - `env/` was abolished and the setting are now under `QligFEP/settings/`;
  - `share/` was renamed to `Qgpu`, to make things more clear.
