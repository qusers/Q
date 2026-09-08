# Two-radius campaign assembly — software validation

The campaign framework writes plans for Q's existing molecular dynamics (MD).
It launches no simulation and submits no high-performance computing (HPC) job.
The unchanged seed matrix and short end-to-end campaign now pass with the existing
SHAKE repair from the modernization branch; see the
[integration evidence](MODERNIZATION_INTEGRATION.md) and
[historical constraint diagnosis](CONSTRAINT_BLOCKER.md).
Do not use software checks as a production-readiness certificate.

## Fixed matrix and profiles

Both profiles require the same sixteen cells: requested radii 10 and 14 angstrom,
charge signs -1 and +1, forward and reverse directions, and replicas 1 and 2 for
each combination. State identities remain charge zero and the declared sign.
Each cell requires a separate, pinned [endpoint-preparation plan](ENDPOINT_PREPARATION.md).

| Profile | Endpoint preparation | Charge ladder | Aggregate duration |
| --- | --- | --- | --- |
| `feasibility_pilot` | 100 picoseconds (ps), including grid seed | Eleven weights 0, 0.1, ..., 1; 30 ps each, with 10 ps designated discard | 6.88 nanoseconds (ns) |
| `software_smoke` | 0.1 ps, including grid seed | Three weights 0, 0.5, 1; 0.1 ps each | 6.4 ps |

Reverse ladders traverse the weights in reverse without swapping state definitions.
The candidate timestep is 1 femtosecond (fs); the 0.5 fs variant preserves physical
durations and doubles step counts. The pilot requires 6,880,000 or 13,760,000
steps respectively. Those are prospective trajectory budgets, not HPC wall times,
allocation approvals or sufficient sampling claims. Failed attempts consume the
experimental allocation and are never silently retried.

All cells use the integrated exterior-Born accounting mode, unchanged real
Lennard–Jones (LJ, repulsion/dispersion) parameters in both charge states, frozen
zero angular offsets and retained restart velocities. Raw/with-Born analysis views
share trajectories; the framework does not duplicate trajectories for a force-free
constant. No ghost endpoints, new angular convention or fitted parameters enter.

## Identity, replication controls and provenance

The checker requires matching isolated builds and the unchanged repository
force-field assets from baseline commit `55aa555c`. Their SHA-256 fingerprints
(SHA-256 is a cryptographic hash algorithm) are pinned in the campaign module.
The probe/library definitions must match across all cells. Radius and orientation
labels must match actual topology headers and preparation commands.

The exact full topology is pinned within every chain. Across independently
initialized replicas, coordinate differences are expected. For this **single
centered probe plus water only**, the checker compares all non-coordinate topology
content, except the header date. It requires the sole probe coordinate to be zero
and all remaining coordinates to be finite, complete water triples. The same-radius
parameter fingerprint, effective radius and frozen offsets must match. This is
not a general permission to ignore coordinate-dependent protein restraints.

Across the matrix, repeated orientation seeds, velocity seeds, prepared coordinate
fingerprints, preparation-report paths or initial restart fingerprints are rejected.
Distinct values are necessary controls against copying a replica, **not proof of
equilibrium independence**. The seed's all-evaluation native diagnostic must pass;
normal termination alone is insufficient. This condition previously blocked
the matrix and is unchanged after the solver repair; no favorable seed was substituted.

## Commands and endpoint transfer

The assembler takes sixteen existing endpoint plans and writes a new directory:

```sh
PYTHONPATH="$PWD/src" python -m QligFEP.charge_campaign assemble \
  /path/to/new-campaign /path/to/endpoint1/plan.json /path/to/endpoint2/plan.json \
  ...remaining-fourteen-endpoint-plan-paths... --profile feasibility_pilot
```

The ellipsis is explanatory, not an executable glob. Supply all sixteen actual
paths. Assembly preserves a failed output directory for diagnosis. Its
`campaign.json` file uses JavaScript Object Notation (JSON), pins endpoint plans,
declares the fixed budget and assigns a unique ladder destination to each cell.
Changing labels, dropping a cell or reducing the reported budget is rejected.

```sh
PYTHONPATH="$PWD/src" python -m QligFEP.charge_campaign inspect /path/to/campaign.json
PYTHONPATH="$PWD/src" python -m QligFEP.charge_campaign stage-ladder \
  /path/to/campaign.json --index 0
```

Inspection is read-only. Staging writes only one ladder and requires its endpoint
preparation to be completely revalidated first; nonexistent future restarts are
not assigned fictitious hashes. The ladder uses schema version 3 with a pinned
endpoint origin, consumes the exact verified final restart, and retains the same
build, topology, charge-state identities, first weight, restraints and velocity
mode. Only step count changes between endpoint preparation and ladder windows.
Every later chain inspection rechecks this origin. Origins pointing to another
ladder are rejected before recursive inspection.

Neither staging nor normal completion sets `equilibrated` or
`independent_replica_established` true. Both are explicitly false. A new output
directory must never be mistaken for evidence of endpoint equilibration.

## Current evidence and remaining work

The software matrix generates sixteen fresh preparations and twenty-step seeds.
After integration of repaired SHAKE, all sixteen cells completed their 0.1 ps
preparation and three 0.1 ps charge windows, with restart provenance and
saved-energy accounting checked. All eleven campaign tests passed with
expected-failure handling disabled before the obsolete markers were removed.
The aggregate 6.4 ps is a software smoke test, not an equilibration experiment.
Short-window analysis correctly reports insufficient sampling, not a usable
confidence interval.

The remaining production work includes multi-observable equilibration/discard/block
sensitivity, between-replica/direction/radius analysis, refreshed same-build hardware
timing and an approved scheduler allocation. No physical correction claim follows
from software matrix consistency alone.
