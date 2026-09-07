# Temporary BAR analysis without exact endpoints

This is an analysis workaround for two-state, linearly mixed Q energy files,
not a change to dynamics, a soft-core implementation, or an analytical boundary
correction. It leaves original files untouched and creates a separate analysis
directory. The normal QFEP header settings, including equilibration discard,
are preserved. It supports full-Hamiltonian (`gas=0`), uncoupled two-state QFEP.

## Why separate this from SCAAS?

Saving unscaled pure-state energies is intentional: QFEP uses the mapping
weights to reconstruct each sampled Hamiltonian. At an exact decoupled endpoint,
a ghost atom can overlap solvent without an interacting-state repulsive force.
Its counterfactual interacting-state energy can consequently become enormous.
That creates an endpoint overlap/numerical problem distinct from whether the
boundary forces, stored boundary energies, and analytical state constants agree.

Exclude complete endpoint ensembles and transitions touching them. Do not clip
energies, delete individual high-energy frames, or replace problematic values
with zero. Rerun QFEP on the retained files; do not merely hide endpoint rows
from an already completed or failed full-ladder calculation.

## Prepare and run

Use a Python environment with this checkout installed (or set `PYTHONPATH` to
its `src` directory). Supply original MD input files containing `[files] energy`
and `[lambdas]`; mapping weights are not inferred from filename digits.

```bash
python -m QligFEP.endpoint_trim prepare /absolute/campaign/qfep.inp \
  --md-input-dir /absolute/campaign/inputfiles \
  --output-dir /absolute/campaign/analysis-interior
```

The output directory must not exist. By default, only exact λ1=0 and λ1=1
windows are excluded: a sampled 0.0001 window is retained. All interior windows
listed in the original QFEP input must be available, uniquely mapped, and
monotonically ordered. At least two interior windows are required. The tool
cannot detect simulations that were never listed in the original input.

For QligFEP's usual layout (QFEP input in `inputfiles/`, energy files in a
replica's `FEP1/298/1/` directory), add `--energy-dir /absolute/path/FEP1/298/1`.
This explicitly selects the working directory used to resolve relative energy
filenames and records it in the manifest; the original QFEP input is unchanged.

For comparisons, explicitly require the same sampled interval in every leg,
replica, and boundary-model variant. For example, add:

```bash
  --lambda-min 0.0001 --lambda-max 0.9999
```

Both interval bounds must actually be sampled; the utility will not silently
substitute more distant windows. These bounds refer to the first mapping weight,
not a ligand identity. Also check that state identities and transformation
directions agree across legs. Changing a previously frozen analysis interval
is a new analysis and should not overwrite or silently replace its manifest.

Run the campaign's validated QFEP executable in the new directory:

```bash
cd /absolute/campaign/analysis-interior
/absolute/path/to/qfep < qfep.inp > qfep.out
python -m QligFEP.endpoint_trim summarize . > summary.json
```

Short symlinks avoid QFEP's filename-length limit. `endpoint-trim.json` records
the source input, MD mappings, retained/excluded files, hashes, and actual
interval/direction. The summarizer checks the prepared input and energy hashes,
BAR row count and λ labels, finite values, and cumulative consistency. Its
result is explicitly marked `full_endpoint_free_energy: false`.

Do not pass this output to the legacy `qligfep_analyze` result reader: its
extraction assumes a final λ1=0 row. Use the dedicated summary. A successfully
parsed result is not proof of overlap, equilibrium sampling, or correct physics;
inspect the QFEP diagnostics and replicate/block stability separately. In
particular, removing exact endpoints does not guarantee sufficient overlap in
the remaining near-endpoint windows.

## Interpretation and analytical constants

The result is F(end mapping) minus F(start mapping), not F(λ1=0) minus
F(λ1=1). Match this interval before subtracting protein and solvent legs.
Cancellation of their missing endpoint contributions must not be assumed.
Do not divide the reported free energy by the retained λ span to extrapolate
to the full transformation. With singular no-softcore interactions, a tiny
omitted λ interval alone does not bound the missing free energy.

For a coordinate-independent state term mixed as

    g(λ) = λ1 g1 + λ2 g2,  with λ1 + λ2 = 1,

its exact contribution over the retained interval is

    [λ2(end) - λ2(start)] (g2 - g1).

The signed prefactor is stored as `state_constant_gap_multiplier`. It is 0.9998
for λ1=0.9999 → 0.0001 and changes sign for the reverse direction. This rule
applies only to genuinely coordinate-independent, linearly mixed terms, not
to the entire free energy or an arbitrary nonlinear correction. Never add a
Born or other state constant twice if it is already included in saved energies.
The rule does not establish that a proposed correction's physics is valid.

## Verification

```bash
python -m pytest -q test/qligfep/test_endpoint_trim.py
QFEP_ENDPOINT_TEST_BINARY=/absolute/path/to/qfep \
  python -m pytest -q test/qligfep/test_endpoint_trim.py
```

The optional native integration uses a two-configuration model with known
partition functions, including deliberately invalid endpoint files that must
never be read. The other tests cover interval selection, direction, preserved
inputs, mapping failures, provenance changes, and malformed/incomplete BAR
output. This is a utility test, not validation of a molecular campaign.
