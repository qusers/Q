"""Prepare a synthetic charge-only probe and run bounded native MD timing checks.

This is not an equilibration, free-energy or production qualification tool.
It calls the existing Qprep/Qdyn programs; it implements no dynamics or sampler.
"""
from __future__ import annotations

import argparse
import json
import math
from pathlib import Path
import platform
import shutil
import subprocess
import time

from . import boundary_native as bn
from .charge_protocol import fingerprint, restart_offsets
from .charge_completion import frames

ROOT = Path(__file__).resolve().parents[2]
PROBE_LIBRARY = """{PRB} !Synthetic steric probe, not a parameterized physical ion
[atoms]
1 C1 CT 0.0
[charge_groups]
C1
"""
PROBE_PDB = "ATOM      1  C1  PRB A   1       0.000   0.000   0.000  1.00  0.00           C\nEND\n"


def _json(path, value):
    with path.open('x') as stream:
        json.dump(value, stream, indent=2, sort_keys=True, allow_nan=False)
        stream.write('\n')


def _execute(binary, directory, args, *, stdin=None):
    start = time.perf_counter()
    # Bound every native invocation; preserve failed inputs/logs for diagnosis.
    with (directory/'native.log').open('x') as log:
        result = subprocess.run([str(binary), *args], cwd=directory, input=stdin,
                                text=True, stdout=log, stderr=subprocess.STDOUT,
                                timeout=60)
    elapsed = time.perf_counter()-start
    if result.returncode:
        raise ValueError(f'Native program failed ({result.returncode}); see {directory}/native.log')
    return elapsed


def prepare(directory, qprep, radius, seed):
    """Generate fresh grid water with Qprep; retain the actual effective radius."""
    if not math.isfinite(radius) or not 8 <= radius <= 22:
        raise ValueError('Bounded probe preparation accepts grid radii 8 through 22 angstrom')
    if abs(radius-round(radius, 3)) > 1e-10:
        raise ValueError('Grid radius must be representable at topology precision (0.001 angstrom)')
    if type(seed) is not int or not 1 <= seed < 100_000_000:
        raise ValueError('Qprep orientation seed must be an integer in [1, 100000000)')
    qprep = qprep.resolve(strict=True)
    qprep_hash = fingerprint(qprep)
    directory.mkdir(parents=True, exist_ok=False)
    for suffix in ('lib', 'prm'):
        shutil.copyfile(ROOT/f'src/QligFEP/FF/AMBER14sb.{suffix}', directory/f'water.{suffix}')
    (directory/'probe.lib').write_text(PROBE_LIBRARY)
    (directory/'probe.pdb').write_text(PROBE_PDB)
    commands = f"""rl water.lib
rl probe.lib
rprm water.prm
rp probe.pdb
set solvent_pack 2.4
set solute_density 0.05794
set random_seed_solute 179857
set random_seed_solvent {seed}
boundary 1 0 0 0 {radius:.6f}
solvate 0 0 0 {radius:.6f} grid HOH
maketop charge_probe
writetop system.top
wp system.pdb y
q
"""
    (directory/'prepare.inp').write_text(commands)
    _execute(qprep, directory, [], stdin=commands)
    log = (directory/'native.log').read_text()
    if 'ERROR' in log.upper() or not (directory/'system.top').is_file():
        raise ValueError(f'Qprep did not produce an error-free topology; see {directory}/native.log')
    lines = (directory/'system.top').read_text().splitlines()
    dimensions = [line for line in lines if '= Total no. of atoms, no. of solute atoms.' in line]
    radii = [line for line in lines if '= Exclusion, solvent radii' in line]
    if len(dimensions) != 1 or len(radii) != 1:
        raise ValueError('Unsupported Qprep topology header')
    atoms, solute = map(int, dimensions[0].split('=')[0].split())
    exclusion, effective = map(float, radii[0].split('=')[0].split())
    if solute != 1 or atoms <= 4 or (atoms-1) % 3 or effective <= 0 or exclusion != radius:
        raise ValueError('Unexpected generated probe dimensions/radii')
    for sign, name in [(-1, 'negative'), (1, 'positive')]:
        (directory/f'{name}.fep').write_text(
            f'[FEP]\nstates 2\n[atoms]\n1 1\n[change_charges]\n1 0 {sign}\n')
    if fingerprint(qprep) != qprep_hash:
        raise ValueError('Qprep binary changed during preparation')
    assets = {p.name: fingerprint(p) for p in sorted(directory.iterdir()) if p.is_file()}
    report = {'schema_version': 1, 'gate': 'probe_preparation_only', 'production_ready': False,
              'grid_radius_angstrom': radius, 'effective_radius_angstrom': effective,
              'exclusion_radius_angstrom': exclusion, 'atoms': atoms, 'waters': (atoms-1)//3,
              'orientation_seed': seed, 'qprep': {'binary': str(qprep), 'sha256': qprep_hash},
              'probe_type': 'existing Amber CT, charge 0 to +/-1; not a physical ion',
              'assets_sha256': assets}
    _json(directory/'prepared.json', report)
    return report


def md_input(radius, steps, timestep, seed, weight):
    """Existing spherical MD; a shared harmonic restraint keeps the probe near zero."""
    return f"""[MD]
steps {steps}
stepsize {timestep}
temperature 298
bath_coupling 10
random_seed {seed}
initial_temperature 298
shake_solvent on
shake_hydrogens on
shake_solute off
lrf off
[cut-offs]
solute_solvent 99
solute_solute 99
solvent_solvent 99
q_atom 99
lrf 99
[sphere]
shell_force 10
shell_radius 0.85
[solvent]
radius {radius}
radial_force 60
polarization on
polarization_force 20
charge_correction on
perstate_polarization on
polarization_adaptation off
perstate_born_correction on
born_dielectric 80
[atom_restraints]
1 0 0 0 10 10 10 0
[intervals]
output 10
non_bond 1
energy 10
[files]
topology system.top
fep charge.fep
final final.re
energy states.en
[lambdas]
{1-weight:.8f} {weight:.8f}
"""


def _check_energies(path, weights, expected_frames):
    """Read only the supported native two-state sequential-record dialect."""
    if sum(1 for _ in frames(path, weights)) != expected_frames:
        raise ValueError('Unexpected number of saved state-energy frames')


def timing(prepared, directory, qdyn, *, sign, weight, steps=1000, timestep=1., seed=112):
    if type(steps) is not int or not 20 <= steps <= 2000 or steps % 10:
        raise ValueError('Timing check requires 20..2000 steps, divisible by 10')
    if sign not in (-1, 1) or weight not in (0., .5, 1.) or timestep not in (.5, 1.):
        raise ValueError('Timing check accepts signs +/-1, weights 0/.5/1, timesteps 0.5/1 fs')
    if type(seed) is not int or not 1 <= seed < 100_000_000:
        raise ValueError('Invalid velocity seed')
    spec_path = prepared/'prepared.json'
    spec = json.loads(spec_path.read_text())
    if spec['gate'] != 'probe_preparation_only' or spec['production_ready'] is not False:
        raise ValueError('Require a probe preparation report')
    for name, checksum in spec['assets_sha256'].items():
        if Path(name).name != name or fingerprint(prepared/name) != checksum:
            raise ValueError('Prepared asset changed or has invalid path')
    qdyn = qdyn.resolve(strict=True)
    engine_hash = fingerprint(qdyn)
    directory.mkdir(parents=True, exist_ok=False)
    shutil.copyfile(prepared/'system.top', directory/'system.top')
    shutil.copyfile(prepared/('positive.fep' if sign == 1 else 'negative.fep'), directory/'charge.fep')
    (directory/'run.inp').write_text(md_input(spec['effective_radius_angstrom'], steps, timestep, seed, weight))
    elapsed = _execute(qdyn, directory, ['run.inp'])
    log = (directory/'native.log').read_text()
    if 'terminated normally.' not in log or 'terminated abnormally' in log:
        raise ValueError('Native run did not terminate normally')
    native = bn.parse(log)
    if (native['meta'][:5] != [spec['atoms'], 1, spec['waters'], 1, 2] or
            native['convention'] != [1] or native['flags'] != [1, 1, 0, 0, 0, 0] or
            native['parameters'][5:7] != [0, 0] or native['water_compatibility'] != [1, 1]):
        raise ValueError('Unexpected native probe model')
    bn._close(native['parameters'][0], spec['effective_radius_angstrom'], 'probe radius')
    for row, q in zip(native['state'], (0, sign)):
        bn._close(row[2], q, 'probe charge')
    if any(row[3] != 0 for row in native['shell']):
        raise ValueError('Timing check requires zero frozen offsets')
    offsets = restart_offsets(directory/'final.re')
    if offsets['atoms'] != spec['atoms'] or offsets['offsets_radians'] != [0]*int(native['meta'][5]):
        raise ValueError('Final dimensions/offsets changed')
    # Q writes steps 10,20,... strictly before the last step.
    frames = (steps-1)//10
    _check_energies(directory/'states.en', [1-weight, weight], frames)
    if fingerprint(qdyn) != engine_hash:
        raise ValueError('Qdyn binary changed during timing check')
    report = {'schema_version': 1, 'gate': 'bounded_md_timing_only', 'production_ready': False,
              'prepared_sha256': fingerprint(spec_path), 'native': native,
              'engine': {'binary': str(qdyn), 'sha256': engine_hash},
              'platform': platform.platform(), 'machine': platform.machine(),
              'steps': steps, 'timestep_fs': timestep, 'sign': sign, 'state2_weight': weight,
              'wall_seconds_including_startup': elapsed, 'saved_frames': frames,
              'ns_per_day_including_startup': steps*timestep/1e6*86400/elapsed,
              'assets_sha256': {p.name: fingerprint(p) for p in sorted(directory.iterdir()) if p.is_file()},
              'limitations': ['unequilibrated grid start; not a production restart',
                              'single short serial run; not HPC performance or convergence',
                              'no proof of binary-to-source provenance',
                              'timing check is not the staged/native campaign restraint gate',
                              'no full-trajectory cutoff/stability audit']}
    _json(directory/'timing.json', report)
    return report


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    commands = parser.add_subparsers(dest='command', required=True)
    prep = commands.add_parser('prepare')
    prep.add_argument('directory', type=Path)
    prep.add_argument('--qprep', type=Path, required=True)
    prep.add_argument('--radius', type=float, required=True)
    prep.add_argument('--seed', type=int, required=True)
    bench = commands.add_parser('timing')
    bench.add_argument('prepared', type=Path)
    bench.add_argument('directory', type=Path)
    bench.add_argument('--qdyn', type=Path, required=True)
    bench.add_argument('--sign', type=int, required=True)
    bench.add_argument('--weight', type=float, required=True)
    bench.add_argument('--steps', type=int, default=1000)
    bench.add_argument('--timestep', type=float, default=1.)
    bench.add_argument('--seed', type=int, default=112)
    args = vars(parser.parse_args())
    command = args.pop('command')
    try:
        report = prepare(**args) if command == 'prepare' else timing(**args)
    except (ValueError, OSError, KeyError, TypeError, subprocess.TimeoutExpired) as error:
        parser.exit(2, f'Probe {command} failed: {error}\n')
    print(json.dumps(report, indent=2, sort_keys=True, allow_nan=False))


if __name__ == '__main__':
    main()
