# Decision needed: angular background-charge convention

Status: evidence-backed proposal, **not implemented or physically validated**.

Post-smooth-pilot update: the [source/response audit](SOURCE_RESPONSE_CONVENTION.md)
quantifies why a single total charge cannot generally replace both the angular
field source and the charging potential in the extended proteins. The compact
enclosed candidate below is not a production recommendation for these geometries.
This is a choice within the charged-boundary investigation, not a proposal for a
new sampler or solvent model. No production target is changed by this document.

Update, 2026-09-09: the [radius/background audit](RADIUS_BACKGROUND_CONVENTION.md)
confirms that the fixed included protein charge is not geometrically enclosed
by the water radius in the completed c-Met/Eg5 campaign. The total-charge
candidate below remains relevant to compact, enclosed systems; it must not be
implemented for these extended proteins merely by substituting the masked
charge sum. The next deliverable is the exterior-response derivation and its
supported geometry, before selecting a replacement production target.

## What the sources and code establish

The [original Surface Constraint All-Atom Solvent (SCAAS) paper](https://doi.org/10.1016/0009-2614(85)87168-2)
defines the angular reference through the dipole response to the solute field
in the complete system (equation 14). Its surface constraint and exterior
continuum contribution have distinct roles. The
[Q methods paper](https://doi.org/10.1016/S1093-3263(99)00012-1)
describes the angular correction as a central-charge approximation for a
nonneutral solute. It also discusses neutralizing distant protein groups and
treating their omitted interactions separately. Neither description establishes
an alchemical atom label as the physical source of the field.

These statements were checked against the supplied local paper transcriptions,
linked in the [initial audit](../../docs/analytical-charge-corrections/INITIAL_AUDIT.md).
The inference is limited: a total-enclosed-charge monopole is a physically
motivated *candidate*, not proof that it provides the correct local field around
an off-center ligand or extended protein. Here monopole means the leading field
contribution determined only by total charge.

In this branch, `wat_shells` derives the angular target from `q_region_state`,
while the Born constant uses `Q_env + q_region_state`. `Q_env` is the nonexcluded,
non-Q solute charge; the Q region is the set designated for state-dependent
treatment. The Q-only angular convention already exists in the parent of
`4ea6c95f` (the endpoint-resolved Hamiltonian commit). It was not introduced by
the recent state-energy bookkeeping changes. Its adequacy in the historical
neutralized-environment workflow does not establish adequacy for retained net
protein charge.

The [native partition check](../../docs/analytical-charge-corrections/PARTITION_AUDIT.md)
shows that relabeling an unchanged +1 atom preserves Born constants but changes
the isolated angular state-energy gap by about 0.896 kilocalories per mole on
one configuration. This is not an estimate of a free-energy error. A small
non-angular numerical mismatch was separately isolated, not hidden.

## Can the frozen scalar offsets represent the missing background?

Write the unclamped target in one shell as

    t_i(q) = theta_base,i - (3/2) A q sin(theta_base,i),
    A = (1 - 1/epsilon) / (rho mu 4 pi r_shell^2).

Here `rho` is water number density, `mu` its molecular dipole magnitude,
`epsilon` the dielectric constant, and `r_shell` the shell reference radius.
The harmonic displacement is `theta_i - t_i(q) + a`, where `a` is the shared
scalar offset. To reproduce the total-charge target with a Q-only target on
the same configuration, a compensating offset would have to satisfy

    a_Q,i - a_total = t_i(q) - t_i(q + Q_env)
                    = (3/2) A Q_env sin(theta_base,i).

The right-hand side varies with rank `i` for a general shell population. One
scalar cannot exactly reproduce it. Clamped targets also require the actual
rankwise target differences; a constant replacement is not generally justified.
The pure native helper regression checks the identity for 101 ranks and
background charges -3, 0 and +3 in an unclamped example. It confirms agreement
at zero background and a nonconstant required shift at nonzero background.
This calculation does not select or fit an offset.

Explicit protein–water forces remain present, but that does not establish that
they make a Q-only restraint equivalent to a total-field restraint. Adaptive
offsets can adjust a mean response; once frozen, they remain ordinary Hamiltonian
parameters. Matching a mean is weaker than matching all target forces and does
not establish the correct charge response. The historical offset discrepancy
still must not be subtracted as an assumed constant free-energy correction.

## Recommended bounded next implementation

Subject to approval, add an **opt-in total-enclosed-charge angular convention**
for controlled validation while retaining the current Q-only default as the
historical control. It would reuse the existing endpoint targets, harmonic
restraint, shell geometry, molecular dynamics (MD), and force/energy bookkeeping.
It must not introduce a new sampler, integrator, thermostat, smoothing function,
radial wall or empirical parameter fitting.

Before trajectories, its native tests must establish:

- Identical results for zero non-Q background charge.
- Angular invariance under fixed-charge Q-region relabeling, with the known
  non-angular control kept separate.
- Consistent included/excluded charge definitions in both angular and Born terms;
  excluded Q atoms cannot be silently included by one and omitted by the other.
- Correct endpoint energy/force mixing and serialization for both charge signs.
- Explicit input, log and manifest identification of the convention, with no
  automatic migration or reinterpretation of historical results/restarts.

The decision is whether to test this minimal candidate, **not** whether it is
already correct. Changing the charge argument changes water forces and sampled
configurations. The current preflight must not accept the new option until its
tests and provenance handling exist. If the candidate needs additional physical
terms or redesigned restraints, pause for a new scope decision.

## What a subsequent physical test must answer

A central, real Lennard–Jones (LJ, repulsion/dispersion) charge-only perturbation
should first test both ligand-charge signs and background-charge signs, with a
zero-background control. Compare radii and matched forward/reverse Hamiltonians
using independently initialized replicas and statistical uncertainty. Test the
chosen model's radius dependence and larger-domain consistency; do not require
positive and negative charging free energies to be equal. Water's molecular
charge geometry and a charged environment need not produce that symmetry, and
the Born environment–perturbation cross term itself is sign dependent.

Then test off-center/distributed environments. Net-charge invariance alone does
not validate those geometries. The exterior monopole approximation, neglected
protein sites, local solvent response, overlap and exact ghost endpoints remain
distinct questions. This option would not remove the existing hard shell-membership
changes or establish exact equilibrium sampling by the existing MD/thermostat.
If those prevent reliable corrections, this candidate must not be declared
validated merely because the charge-partition checks pass. No charged result,
binding affinity or experimental value
may be used to tune the candidate. A measured compute budget and user approval
are required before substantial high-performance computing (HPC) calculations.

If this target change is not approved, the existing convention remains intact.
It can still be studied as a defined model, but the current evidence does not
support declaring it a generally valid correction for retained charged proteins.

Verification checkpoint: the expanded native helper assertions pass; the full
focused suite remains **107 passed, 2 skipped, 2 expected failures in 14.72 seconds**.
Only tests and documentation changed in this decision checkpoint. The existing
engine target, historical data and experiment inputs were not modified, and no
HPC job was submitted. `git diff --check` passed.
