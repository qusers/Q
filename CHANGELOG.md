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
- Three CLIs implemented: `qligfep`, `setupFEP`, and `qmapfep`.
  - TODO: `setupFEP` is currently setup in a way that it creates only the **protein** and the **water** legs. In the future, it would be good to have it flexible so that it can also create the **vacuum** legs as well.
- q6's makefile no longer moves the executables to Q/ and instead keeps then within `src/q6`. Useful for future improvements where the user will get everything necessary to get started by doing pip install.
- openff2Q now implemented within a CLI as well. Still need to update the module for supporting the previous ForceFields as well (e.g.: charmm, opls...).
- openff2Q calculates atom masses directly from the package instead of the previous hardcoded dictionary.
- created the `pdb_to_amber.py` module, accessible through the CLI as `pdb2amber` to programatically rename the aminoacids according to the force-field conventions.
  - TODO: in the script there's also a bit for identifiying cystein bonds. Should move this somewhere else in the future.
- `analyze_FEP.py` now accessible through the CLI with the command `qligfep_analyze`.
- Switched the SHAKE algorithm to `off` on [equilibrations 1 to 4](https://github.com/qusers/Q/commit/4758c2b010cead0d40686b7c608405a4c576ab05). TODO: benchmark SHAKE behavior when ON from equilibration 4 as well.
- Changed the timelimit for jobs on Tetralith from 24:00:00 to 12:00:00 to allow our jobs to be used as ["backfill" jobs](https://www.nsc.liu.se/support/batch-jobs/tetralith/adapt-to-workflow/#:~:text=If%20possible%2C%20structure%20your%20jobs%20so%20that%20the%20time%20limit%20is%20between%2030%20minutes%20and%2024%20hours.%20Such%20jobs%20can%20often%20be%20started%20using%20%22backfill%22%20(using%20nodes%20that%20would%20otherwise%20be%20idle).).



## TODO:
- [x] add a MANIFEST.in file to include non-python files in the package.
- [x] Update Qgpu with the correct imports as well.
- [  ] Implement Softcore to the package. ⚠️ Would take a long time ⚠️
  - Can take a look at the implementation of this in other packages and translate it to fortran code
- [x] Make a script detecting the cystein bonds so that the user doesn't need to [hardcode](https://github.com/GPCR-ModSim/qligfep-benchmark/blob/main/inputfilegen/OPLS2015/Thrombin/setup.py) them anymore.
- [  ] Numpy seems to have a support for compiling fortran code. Could think about making a wrapper so that thes functions can be directly called with python. (for now, lower priority...)
- [x] Setup a proper logging to the package; for now this was (roughly) done with [loguru](https://github.com/Delgan/loguru).
- [x] Make current `except` statements in the code more specific. When not specific, these statements can also catch, for example, `KeyboardInterrupt`, triggered by a `Ctrl + c` in the terminal. (but ongoing...)
- [  ] Currently, QmapFEP makes a `self.ligands` dictionary with the charges as keys. However, that's problematic because the charges were hard-coded... I added the charges to the MoleculePool but the code doesn't use that yet.
- [  ] Add test modules for the current refactored code.

## TODO:
- qpenff2Q module has a structure that would be useful to be made into a base class. E.g.: we could use it as a child class for handling routine operations with `.sdf` files.
- QmapFEP not working when .sdf file contains more than 1 unique formal charges.
