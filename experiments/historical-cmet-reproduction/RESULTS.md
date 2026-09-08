# Historical CMET interior BAR reproduced

Date: 2026-09-07. **Numerical reproduction gate passed.** This is not physical
validation of the boundary model or a full-endpoint free-energy calculation.

## Scope and result

The [predeclared selection](selection.json) is the first edge in the archived
CMET campaign list, `CHEMBL3402741_400 -> CHEMBL3402744_300`, replica 1, both
directions and both legs. No selection based on experiment, overlap, or a
favorable result was performed. No new molecular simulation was run.

| Reduction | Reproduced interior BAR | Archived cumulative difference | Maximum error across all 99 rows |
| --- | ---: | ---: | ---: |
| Forward water | 70.944 | 70.944 | 0.001 |
| Forward protein | 34.592 | 34.592 | 0.001 |
| Reverse water | 128.287 | 128.288 | 0.001 |
| Reverse protein | 162.630 | 162.630 | 0.001 |

All entries are kcal/mol. All four reductions pass the fixed 0.002 kcal/mol
tolerance at every retained row (396 row comparisons, 392 adjacent pairs).
The reference values are differences between the archived full-ladder cumulative
BAR columns at lambda1=0.999 and 0.001, not an independently published trimmed
reference. The new calculations remove both exact-endpoint files and rerun
the retained ladder from zero origin.

Preserved settings:

- lambda1=0.999 -> 0.001: 99 retained windows, 98 pairs per reduction;
- 499 saved frames per window, first 100 discarded: 399 analyzed;
- QFEP kT=0.592 kcal/mol, not recomputed from the directory name `298`;
- QFEP state-2 alpha=100 kcal/mol and all other input header values unchanged;
- the historical safeguarded BAR source, rebuilt locally rather than modifying
  or executing the archived Linux binary.

The alpha constant contributes **99.8 kcal/mol to each listed leg** and cancels
in protein-minus-water differences. Thus those large individual numbers must
not be mistaken for analysis-offset-free physical charging free energies.
It is preserved here for faithful numerical reproduction, not adopted as an
additional physical correction.

## Born convention and interpretation

The campaign manifest identifies dynamics engine `55aa555c`, endpoint-resolved
polarization enabled and integrated Born disabled. Its archived analysis uses
`ke=332.0637`, epsilon=80, and the rounded charges/radii in native logs. That
convention was independently reconstructed; it is deliberately **not replaced
retroactively** by the new topology-constant fix when reproducing old results.

The logged environment charges are 0 (water) and +3 (protein); radii are 24.730
and 23.300 A. Forward Q charges are +1.000 and -0.004, exchanged in reverse.
The full state-gap Born correction is +42.796318 kcal/mol protein-minus-water
forward. Multiplying by the retained span 0.998 gives +42.710726 kcal/mol.

| Direction | Raw protein minus water | Historical Born over retained interval | Corrected interior difference |
| --- | ---: | ---: | ---: |
| Forward | -36.352 | +42.710726 | +6.358726 |
| Reverse | +34.343 | -42.710726 | -8.367726 |

These are single-replica, truncated quantities with no uncertainty estimate or
experimental comparison. The forward/reverse sum is -2.009 kcal/mol, but
**do not label it same-Hamiltonian closure**: decoding the third record of each
`eq5.re` shows different frozen polarization offsets between directions for
both legs. Their exact single-precision values are retained in [result.json](result.json).
No claim is made here about how much of the directional difference those
offsets explain.

The actual FEP inputs and setup metadata declare `softcore_method gapsys`,
with state-specific softcore assignments. This is archival provenance, not
evidence that the present user's no-softcore system is identical. A method name
alone does not establish effective softening of interacting-state terms.

Successful reproduction does not establish equilibrium sampling, statistical
independence, BAR overlap, the canonical validity of ranked SCAAS, or cancellation
of omitted endpoint caps. None of those gates is claimed to pass here.

## Integrity and repeatability

All **872 downloaded files (63,998,149 bytes)** match an aggregate SHA-256
inventory calculated read-only on Snellius. This includes mappings, energy files,
archived outputs, boundary logs, restarts, topologies, and analysis provenance.
The archived QFEP binary MD5 also matches its historical validation record.
Every stored state mapping in every retained frame was checked against the
original MD inputs, and every retained serialized energy component was finite.
The source and binary hashes, constants, offsets, and row-error summaries are
recorded in [result.json](result.json).

Raw data, the archived binary, local build products, and full rerun outputs are
ignored by Git. They remain under `raw/`, `engine-build/`, and `analysis/run1/`;
no source data were overwritten. The development worktree's QFEP is unchanged.
The preparation utility gained a tested `--energy-dir` option for QligFEP's
separate `inputfiles/` and replica-energy directories, avoiding rewritten inputs.

From this case directory, with this worktree on `PYTHONPATH`:

```bash
# Only on a fresh checkout without raw/; requires read access via ssh mysnellius.
python fetch_reference.py

mkdir -p engine-build
cd engine-build
gfortran-11 -O3 -cpp -std=legacy -ffree-line-length-none -DG95=1 \
  ../raw/engine/sizes.f90 ../raw/engine/nrgy.f90 ../raw/engine/misc.f90 \
  ../raw/engine/parse.f90 ../raw/engine/mpiglob.f90 ../raw/engine/qfep.f90 -o qfep
cd ..
python reproduce.py --output-dir analysis/run2
```

Always choose a new output directory. The fetch and reproduction scripts refuse
to overwrite previous reference/results directories. Failed attempts, if any,
should be preserved separately; a timeout is not a passing result.

The regression entry point is `test/qligfep/test_historical_cmet_reproduction.py`.
Its native reference test is opt-in via `Q_RUN_HISTORICAL_CMET=1`; ordinary
parser and record-integrity tests do not need the raw data or a network connection.

Final verification: **52 tests passed, none skipped**, with the historical native
integration and native Qdyn/QFEP tests enabled. `git diff --check` also passed.

## Next development gate

The numerical baseline is now reproducible. Next, define and test the smallest
repair of ranked-SCAAS force/energy consistency at shell crossings, with fixed,
explicitly matched boundary parameters. Reuse the existing discontinuity audits;
do not repeat the neutral-droplet campaign or fit a new correction to these
single-edge values. Full-endpoint ghost overlap remains a separate issue.
