# Frozen offsets are Hamiltonian parameters

In Q's existing molecular dynamics (MD), each Surface Constraint All-Atom Solvent
(SCAAS) shell has an angular offset `theta_corr`. Freezing adaptation does not
remove the offset: it fixes a parameter in the potential and its forces. Different
frozen offsets generally mean different Hamiltonians, even with identical charge
states, topology and boundary radius.

## Rechecked historical evidence

The [preserved historical reduction](../../experiments/historical-cmet-reproduction/result.json)
contains one replicate of a protein/water perturbation in each direction. The
archived raw data were rechecked read-only against the pinned inventory: 872 files,
63,998,149 bytes, inventory SHA-256 (a cryptographic content fingerprint)
`e4aaaa67808158cc7101bf87277d5b834fdd8f7189c2d4282c27357c6dd13e91`.
The four `eq5.re` restart offset records exactly match the recorded reductions.

The following differences are reverse minus forward, with shell 1 outermost:

| Leg | Shell 1, radians (degrees) | Shell 2, radians (degrees) | Shell 3, radians (degrees) |
| --- | ---: | ---: | ---: |
| Water | +0.01314342 (+0.75306°) | -0.01276624 (-0.73145°) | -0.00095189 (-0.05454°) |
| Protein | +0.01161015 (+0.66521°) | +0.00334537 (+0.19168°) | -0.00370777 (-0.21244°) |

The historical interior-only, Born-corrected protein-minus-water results are
+6.358726 and -8.367726 kilocalories per mole in the two directions. Their sum,
-2.009 kilocalories per mole, is **not a same-Hamiltonian closure residual**.
The offset mismatch is a confound; its contribution to that number has not been
measured. One replicate, without an uncertainty/convergence assessment, also
does not establish statistical significance. The historical inputs configure
Gapsys softcore; they must not be conflated with the separate no-softcore ghost
overlap problem described for other calculations.

## Why there is no simple offset subtraction from a free energy

For one shell with force constant `k`, offset `a`, observed angles `theta_i` and
pure-state rank targets `t_s,i`, the implemented harmonic restraint is

    U_s(a) = (k/2) sum_i (theta_i - t_s,i + a)^2.

On the same configuration, changing the offset by `delta` gives exactly

    U_s(a+delta) - U_s(a)
      = k delta sum_i(theta_i - t_s,i + a) + (k N/2) delta^2.

Consequently the change in the state-2 minus state-1 energy gap is

    [U_2-U_1](a+delta) - [U_2-U_1](a)
      = k delta sum_i(t_1,i - t_2,i).

These identities also hold when the target angles are clamped. The pure native
[boundary-function tests](../../test/q6/test_boundary_corrections.f90) check the
energy, gap and angular-gradient identities across observed angles and clamped
and unclamped targets. They do not measure a historical free-energy shift.

The gap identity eliminates the observed angles for a fixed shell population,
but Q's ranked target list depends on the instantaneous number `N` of waters in
that shell. More importantly, the single-state energy shift depends on the
configuration and changes the ensemble. It is therefore not generally a
coordinate-independent correction like the fixed-radius Born term. Subtracting
an average gap shift from the historical estimate is not justified by this
identity. A retrospective common-Hamiltonian estimate would require adequate
configuration data, cross-evaluation and validated reweighting overlap; the
saved state-energy records alone do not establish that.

## Minimal native restart safety change

Previously, a missing offset record or mismatched shell count produced a warning
and reset every offset to zero, even with frozen per-state polarization. The
[native loader](../../src/q6/md.f90) now stops in that frozen mode instead of
silently changing the Hamiltonian. Incomplete or nonfinite offset payloads also
stop with an explicit diagnostic. A fresh run without a restart still initializes
zero offsets, as before. The adaptive legacy fallback for missing records or
mismatched counts remains unchanged.

The [native restart tests](../../test/q6/test_frozen_offset_restart.py) verify
bit-for-bit preservation of valid nonzero offsets in a short MD restart, rejection
before energy sampling for five invalid-record cases, and both adaptive fallback
cases. This does not change a force law, integrator, thermostat, sampler, or
boundary target for a valid input. These tests ran with the serial GNU Fortran 11
build; other compiler and distributed-execution builds remain to be verified.

## Required production control

For each physical system/radius, declare one frozen offset vector and keep its
binary representation identical across all lambda windows, directions, charge
signs and independent replicas intended to estimate the same Hamiltonian.
Independent coordinates and velocities do not require independently chosen
Hamiltonian parameters. Water and protein legs may use distinct declared vectors
because they are distinct systems; a radius comparison must explicitly declare
its offset-selection rule, not assume shell offsets are transferable.

- Record the common source restart and offset values; check the input and final
  restart records for every production window. Require adaptation to be off.
- Check the full boundary definition too: equal shell counts alone do not prove
  equal radii, shell widths, force constants, charge conventions or atom mapping.
- Do not choose offsets to improve charged, binding or experimental agreement.
  Separately adapted forward/reverse targets cannot be called a convergence test.
- Treat old interior-only Bennett acceptance ratio (BAR) estimates as truncated
  intervals. For new charge-only validation keep identical real Lennard–Jones
  (LJ, repulsion/dispersion) interactions in both states so no disappearing-atom
  endpoint is introduced.

The native guard enforces only restart validity, not cross-run parameter matching.
An experiment-level manifest/checker is still required before high-performance
computing (HPC) production. No historical result has been retuned, and no HPC
job was submitted for this audit.

Verification checkpoint: **72 passed, 2 skipped, 2 expected failures in 9.25
seconds** across the focused native and analysis suite, after rebuilding Qdyn
with `make -C src/q6 qdyn FC=gfortran-11`. The eight new native restart tests all
ran. The optional analysis-runtime skips and the explicitly retained non-angular
partition failures are unchanged from the [partition audit](PARTITION_AUDIT.md).
`git diff --check` passed.
