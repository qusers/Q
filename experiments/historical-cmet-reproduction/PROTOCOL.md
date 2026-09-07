# Existing-data CMET reproduction (2026-09-07)

Select the first edge in the archived CMET campaign list:
`CHEMBL3402741_400 -> CHEMBL3402744_300`, and its reverse. Use replica 1,
protein and water: four reductions. Selection is by campaign ordering, not
agreement with experiment, overlap, or closure. No experimental comparison or
new molecular dynamics is part of this gate.

Historical campaign on `mysnellius`:
`/projects/prjs2157/charge-change/experiments/cmet-full-pspol-v1`.
Its manifest identifies the Qdyn engine as `55aa555c`, endpoint-resolved
polarization enabled and integrated Born disabled (Born applied post hoc).

Download only the selected mappings, state/energy files, archived QFEP outputs,
relevant topology/preparation provenance, and the archived safeguarded-QFEP
source/validation record. Remote operations must be read-only. Preserve the
downloaded originals separately from local build and analysis products.

Rebuild that archived QFEP source locally; do not replace the development
engine's QFEP or the archived binary. This reproduces a solver implementation,
not a bit-identical Linux executable. Record source and binary hashes.

Run the new endpoint-trim preparation on each reduction, preserving the
original QFEP analysis parameters and sampled MD lambda mappings. Use the
historical comparison interval lambda1=0.999 to 0.001 if both bounds are sampled;
otherwise stop and document the actual ladder rather than silently change it.
Retain every intervening window and all frames after the original discard.

Reference: extract the difference of cumulative BAR values at those two lambda
rows from each archived safeguarded full-ladder QFEP output. Label this explicitly
as a difference of archived printed cumulative values, not a newly discovered
standalone trimmed reference. Compare every available retained cumulative row
after subtracting the archived 0.999 origin, not only the final total.

Gate: absolute agreement within 0.002 kcal/mol for each row/reduction (allowing
printed-column rounding and the historical solver's documented tolerance).
Reject missing/incomplete/nonfinite outputs and mapping mismatches. Do not clip
energies, skip interior pairs, or interpret solver timeouts as success.

Report raw protein-minus-water reductions and independently reproduce the
archived post-hoc Born convention over the same interval. Separate the physical
claim from reproduction: a match does not validate continuum physics, canonical
sampling, overlap, or convergence. A single replica cannot estimate uncertainty
or establish a forward/reverse closure problem. Check whether the forward and
reverse frozen polarization offsets actually define the same Hamiltonians
before interpreting their sum as thermodynamic closure.

If the new utility cannot consume the archived layout, implement the smallest
explicit, tested compatibility change; never alter the downloaded input silently.
