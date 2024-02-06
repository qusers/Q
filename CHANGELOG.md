# CHANGELOG - QligFEP refactoring:
## New packaging structure:
```text
Q/
├── LICENSE
├── pyproject.toml
├── setup.py # needed for compiling Q
├── README.md
├── src/
│   ├── q6/
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
- Changed a statements like `== True` to `is True`. Also fixed statements such as `!= None` and `type(<object>)` == list to `isinstance()`.
- Command line interface (CLI) moved to `QligFEP/CLI`.
- All `sys.path.insert` statements removed and substituted by relative imports within the package.
- Moved the IO functions `pdb_parse_out`, `pdb_parse_in` to the module `pdb_utils` to avoid circular imports within the package.
- Three CLIs currently implemented: `qligfep`, `setupFEP`, and `qmapfep`. `TODO: will probably make it capitalized in the future.`
- q6's makefile no longer moves the executables to Q/ and instead keeps then within `src/q6`. Useful for future improvements where the user will get everything necessary to get started by doing pip install.
- openff2Q now implemented within a CLI as well. Still need to update the module for supporting the previous ForceFields as well.



## TODO:
- [  ] add a MANIFEST.in file to include non-python files in the package.
- [  ] Update Qgpu with the correct imports as well.
- [  ] Unify both IO modules from Qgpu and QligFEP to a single one that can interact with both implementations [`q6` and `Qgpu` backends].
- [  ] Get Lomap running and compare it to QmapFEP.
- [  ] Implement Softcore to the package. ⚠️ Might take a long time ⚠️
  - Can take a look at the implementation of this in other packages and translate it to fortran code
  - Testing: Would run simulations with & without softcore ➡️ How do values compare? 
- [  ] Make a script detecting the cystein bonds so that the user doesn't need to [hardcode](https://github.com/GPCR-ModSim/qligfep-benchmark/blob/main/inputfilegen/OPLS2015/Thrombin/setup.py) them anymore
- [  ] Numpy seems to have a support for compiling fortran code. Could think about making a wrapper so that thes functions can be directly called with python. (for now, lower priority...)
- [x] Setup a proper logging to the package. An easy option would be to use [loguru](https://github.com/Delgan/loguru), which I personally find very cool 🙂.
- [  ] Make current `except` statements in the code more specific. When not specific, these statements can also catch, for example, `KeyboardInterrupt`, triggered by a `Ctrl + c` in the terminal.
- [  ] Currently, QmapFEP makes a `self.ligands` dictionary with the charges as keys. However, that's problematic because the charges were hard-coded... I added the charges to the MoleculePool but the code doesn't use that yet.

## Questions:
- QmapFEP not working when .sdf file contains more than 1 unique formal charges.
