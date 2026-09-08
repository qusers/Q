"""Validate a planned restart chain and launch at most one realized native window.

No scheduler submission, automatic failed-job retry, or production qualification.
Existing attempt/output files are preserved. Every successor needs a revalidated
completion receipt for its predecessor; future restart hashes are never invented.
"""
from __future__ import annotations

import argparse
from datetime import datetime, timezone
from decimal import Decimal
import json
import math
import os
from pathlib import Path
import platform
import subprocess
import time

from . import charge_protocol as cp, charge_completion as completion
from . import charge_build

RESERVED = ('charge-started.json', 'charge-preflight.json', 'charge-completed.json',
            'charge-failed.json', 'charge-native.log')


def _write(path, value):
    with path.open('x') as stream:
        json.dump(value, stream, indent=2, sort_keys=True, allow_nan=False)
        stream.write('\n')


def inspect_plan(path):
    path = path.resolve()
    spec = json.loads(path.read_text())
    keys = {'schema_version', 'engine', 'series', 'initial_restart', 'initial_restart_sha256',
            'build_report', 'build_report_sha256'}
    endpoint = spec.get('schema_version') == 2 and spec.get('purpose') == 'endpoint_preparation'
    if endpoint:
        keys |= {'purpose', 'preparation'}
    if set(spec) != keys or (not endpoint and spec['schema_version'] != 1):
        raise ValueError('Unsupported chain plan schema')
    build_path = (path.parent/spec['build_report']).resolve()
    if cp.fingerprint(build_path) != spec['build_report_sha256']:
        raise ValueError('Build report hash mismatch')
    build = charge_build.validate(build_path)
    engine = spec['engine']
    if set(engine) != {'binary', 'sha256', 'source_commit'}:
        raise ValueError('Require pinned engine binary and declared source commit')
    binary = (path.parent/engine['binary']).resolve()
    if cp.fingerprint(binary) != engine['sha256']:
        raise ValueError('Engine hash mismatch')
    if len(engine['source_commit']) != 40 or any(c not in '0123456789abcdef' for c in engine['source_commit']):
        raise ValueError('Require full source commit declaration')
    if (binary != (build_path.parent/build['binaries']['qdyn']['path']).resolve() or
            engine['sha256'] != build['binaries']['qdyn']['sha256'] or engine['source_commit'] != build['source_commit']):
        raise ValueError('Engine declaration does not match isolated build record')
    initial_path = (path.parent/spec['initial_restart']).resolve()
    if cp.fingerprint(initial_path) != spec['initial_restart_sha256']:
        raise ValueError('Initial restart hash mismatch')
    initial = cp.restart_offsets(initial_path)
    series = spec['series']
    if set(series) != {'id', 'system', 'sign', 'direction', 'replica', 'born_mode', 'apply_born_posthoc', 'windows'}:
        raise ValueError('Unsupported chain series declaration')
    if not series['id'] or not series['system'] or type(series['replica']) is not int or series['replica'] < 1:
        raise ValueError('Require series/system identity and positive replica number')
    mode = series['born_mode']
    if mode not in ('integrated', 'posthoc', 'control') or type(series['apply_born_posthoc']) is not bool or series['apply_born_posthoc'] != (mode == 'posthoc'):
        raise ValueError('Invalid or double-counted Born accounting mode')
    if type(series['sign']) is not int or series['sign'] not in (-1, 1):
        raise ValueError('Require charge sign +/-1')
    if series['direction'] not in ('forward', 'reverse') or len(series['windows']) < (1 if endpoint else 2):
        raise ValueError('Require a forward/reverse chain including both endpoints')
    protected = {path, binary, initial_path, build_path}
    protected.update((build_path.parent/'native-source.tar', build_path.parent/'build.log'))
    protected.update((build_path.parent/name).resolve() for name in build['source_files_sha256'])
    protected.update((build_path.parent/item['path']).resolve() for item in build['binaries'].values())
    outputs, directories, windows = set(), set(), []
    predecessor = initial_path
    for item in series['windows']:
        if set(item) != {'input', 'sha256', 'assets_sha256'} or set(item['assets_sha256']) != {'topology', 'fep'}:
            raise ValueError('Pin each input and topology/FEP; future restart hashes are not input assets yet')
        inp = (path.parent/item['input']).resolve()
        if cp.fingerprint(inp) != item['sha256']:
            raise ValueError('Planned input hash mismatch')
        window = cp.inspect_definition(inp, mode, initial['atoms'])
        if window['assets_sha256'] != item['assets_sha256']:
            raise ValueError('Planned topology/FEP hash mismatch')
        if Path(window['paths']['restart']) != predecessor:
            raise ValueError('Restart must be initial asset or immediately preceding final output')
        if inp.parent in directories:
            raise ValueError('Each chain window needs its own launch directory')
        directories.add(inp.parent)
        protected.update((inp, Path(window['paths']['topology']), Path(window['paths']['fep'])))
        targets = {Path(window['paths'][key]) for key in ('energy', 'final')}
        if len(targets) != 2 or any(p.parent != inp.parent for p in targets):
            raise ValueError('Distinct energy/final outputs must stay in their window directory')
        reserved = {inp.parent/name for name in RESERVED}
        if targets & reserved or outputs & (targets | reserved):
            raise ValueError('Chain output/receipt collision')
        outputs.update(targets | reserved)
        predecessor = Path(window['paths']['final'])
        windows.append(window)
    if outputs & protected:
        raise ValueError('Chain output would overwrite a protected input or receipt name')
    first = windows[0]
    for window in windows:
        signature = json.loads(json.dumps(window['signature']))
        if endpoint:
            signature['md']['steps'] = first['signature']['md']['steps']
        if (signature != first['signature'] or window['states'] != first['states'] or
                window['assets_sha256'] != first['assets_sha256']):
            raise ValueError('Within-chain Hamiltonian, topology/FEP or simulation settings differ')
        if list(map(Decimal, window['q_region_charges'])) != [0, series['sign']]:
            raise ValueError('Require charge-only 0 to declared sign')
    weights = [Decimal(w['lambdas'][1]) for w in windows]
    ends = (0, 1) if series['direction'] == 'forward' else (1, 0)
    orientation = 1 if series['direction'] == 'forward' else -1
    if endpoint:
        if any(w != ends[0] for w in weights) or first['signature']['velocity_initialization'] != 'restart':
            raise ValueError('Endpoint preparation must retain restart velocities at its fixed starting endpoint')
        from .charge_endpoint import validate_origin
        origin = validate_origin(path, spec, build_path, build, windows)
    elif (weights[0], weights[-1]) != ends or any(orientation*(b-a) <= 0 for a, b in zip(weights, weights[1:])):
        raise ValueError('Require complete monotonic charge-only chain')
    return {'schema_version': spec['schema_version'], 'gate': 'planned_chain_consistency_passed', 'production_ready': False,
            'purpose': 'endpoint_preparation' if endpoint else 'charge_ladder',
            'preparation_origin': origin if endpoint else None,
            'native_trace_required': 'subroutine write_charge_trace' in (build_path.parent/'src/q6/md.f90').read_text(),
            'plan_path': str(path), 'plan_sha256': cp.fingerprint(path),
            'engine': {**engine, 'binary': str(binary)}, 'series': {**series, 'windows': windows},
            'build_report': str(build_path), 'build_report_sha256': spec['build_report_sha256'],
            'driver_files_sha256': {name: cp.fingerprint(Path(__file__).with_name(name)) for name in
                                   ('charge_chain.py', 'charge_build.py', 'charge_protocol.py',
                                    'charge_completion.py', 'boundary_native.py', 'endpoint_trim.py',
                                    'charge_endpoint.py', 'charge_probe.py', 'charge_diagnostics.py')},
            'initial_restart': str(initial_path), 'initial_restart_sha256': spec['initial_restart_sha256'],
            'initial_offsets': initial, 'total_steps': sum(int(w['signature']['md']['steps']) for w in windows),
            'limitations': ['future restart contents are not validated until realized',
                            'single-chain consistency does not prove campaign completeness or independence',
                            'recorded source/build/launch evidence is not a signed attestation',
                            'no equilibration, statistical or physical qualification']}


def _realized(plan, index, predecessor_hash):
    definition = plan['series']['windows'][index]
    window = cp.inspect_window(Path(definition['input']), plan['series']['born_mode'])
    declared = {k: v for k, v in window.items() if k != 'restart'}
    declared['assets_sha256'] = {k: v for k, v in window['assets_sha256'].items() if k != 'restart'}
    if declared != definition or window['assets_sha256']['restart'] != predecessor_hash:
        raise ValueError('Realized window differs from plan or verified predecessor')
    if window['restart']['offset_record_sha256'] != plan['initial_offsets']['offset_record_sha256']:
        raise ValueError('Realized restart changed frozen offsets')
    return window


def progress(plan):
    """Recheck completed outputs, never infer that an incomplete attempt is dead."""
    predecessor_hash = plan['initial_restart_sha256']
    completed = 0
    pending = False
    for index, definition in enumerate(plan['series']['windows']):
        directory = Path(definition['input']).parent
        receipt_path = directory/'charge-completed.json'
        if not receipt_path.exists():
            future_artifacts = [directory/name for name in RESERVED]+[
                Path(definition['paths'][key]) for key in ('energy', 'final')]
            if pending and any(p.exists() or p.is_symlink() for p in future_artifacts):
                raise ValueError('Attempt exists beyond an incomplete predecessor')
            pending = True
            continue
        if pending:
            raise ValueError('Completion receipt exists beyond an incomplete predecessor')
        receipt = json.loads(receipt_path.read_text())
        if receipt['plan_sha256'] != plan['plan_sha256'] or receipt['window_index'] != index:
            raise ValueError('Completion receipt belongs to a different plan/window')
        start_path, preflight_path = directory/'charge-started.json', directory/'charge-preflight.json'
        if (receipt['started_sha256'] != cp.fingerprint(start_path) or
                receipt['preflight_sha256'] != cp.fingerprint(preflight_path)):
            raise ValueError('Launch/preflight receipt changed')
        window = _realized(plan, index, predecessor_hash)
        preflight = json.loads(preflight_path.read_text())
        if preflight != {'plan_sha256': plan['plan_sha256'], 'window_index': index, 'window': window}:
            raise ValueError('Realized preflight does not match actual window')
        started = json.loads(start_path.read_text())
        if (started['plan_sha256'] != plan['plan_sha256'] or started['window_index'] != index or
                started['engine'] != plan['engine'] or started['cwd'] != str(directory) or
                started['driver_files_sha256'] != plan['driver_files_sha256'] or
                started['command'] != [plan['engine']['binary'], Path(window['input']).name]):
            raise ValueError('Launch declaration differs from planned invocation')
        result = completion.validate_window(window, plan['series']['born_mode'], directory/'charge-native.log')
        if plan['native_trace_required'] and result['trajectory_diagnostics'] is None:
            raise ValueError('Retained native source requires trajectory diagnostics; trace is missing')
        if result != receipt['result']:
            raise ValueError('Completed native outputs changed since their receipt')
        if (directory/'charge-failed.json').exists():
            raise ValueError('Conflicting failure and completion receipts')
        predecessor_hash = result['final_sha256']
        completed += 1
    return completed, predecessor_hash


def run_next(path, *, max_steps=2000, timeout=60.):
    if type(max_steps) is not int or max_steps < 1 or not math.isfinite(timeout) or timeout <= 0:
        raise ValueError('Declare positive per-invocation step and time budgets')
    plan = inspect_plan(path)
    index, predecessor_hash = progress(plan)
    if index == len(plan['series']['windows']):
        return {'gate': 'chain_outputs_consistent', 'production_ready': False, 'completed_windows': index}
    window = _realized(plan, index, predecessor_hash)
    directory = Path(window['input']).parent
    if int(window['signature']['md']['steps']) > max_steps:
        raise ValueError('Next window exceeds this invocation step budget')
    protected_outputs = [directory/name for name in RESERVED]+[Path(window['paths'][key]) for key in ('energy', 'final')]
    if any(p.exists() or p.is_symlink() for p in protected_outputs):
        raise ValueError('Existing attempt/output: inspect process or scheduler; never automatically overwrite or retry')
    command = [plan['engine']['binary'], Path(window['input']).name]
    started = {'plan_sha256': plan['plan_sha256'], 'window_index': index, 'command': command,
               'cwd': str(directory), 'engine': plan['engine'], 'wrapper_pid': os.getpid(),
               'driver_files_sha256': plan['driver_files_sha256'],
               'hostname': platform.node(), 'started_utc': datetime.now(timezone.utc).isoformat(),
               'scheduler_context': {key: os.environ[key] for key in
                                     ('SLURM_JOB_ID', 'SLURM_ARRAY_JOB_ID', 'SLURM_ARRAY_TASK_ID') if key in os.environ},
               'max_steps': max_steps, 'timeout_seconds': timeout}
    # Exclusive creation is the launch claim. It is retained on any failure,
    # including a crash before a child PID or completion receipt is written.
    _write(directory/'charge-started.json', started)
    start_time = time.perf_counter()
    try:
        _write(directory/'charge-preflight.json', {'plan_sha256': plan['plan_sha256'], 'window_index': index, 'window': window})
        with (directory/'charge-native.log').open('x') as stream:
            result = subprocess.run(command, cwd=directory, stdout=stream, stderr=subprocess.STDOUT, timeout=timeout)
        if result.returncode != 0:
            raise ValueError(f'Native child exited with status {result.returncode}')
        if inspect_plan(path) != plan:
            raise ValueError('Plan, input or engine changed during execution')
        output = completion.validate_window(window, plan['series']['born_mode'], directory/'charge-native.log')
        if plan['native_trace_required'] and output['trajectory_diagnostics'] is None:
            raise ValueError('Retained native source requires trajectory diagnostics; trace is missing')
        receipt = {'plan_sha256': plan['plan_sha256'], 'window_index': index,
                   'started_sha256': cp.fingerprint(directory/'charge-started.json'),
                   'preflight_sha256': cp.fingerprint(directory/'charge-preflight.json'),
                   'wall_seconds': time.perf_counter()-start_time, 'result': output}
        _write(directory/'charge-completed.json', receipt)
    except Exception as error:
        _write(directory/'charge-failed.json', {'plan_sha256': plan['plan_sha256'], 'window_index': index,
                                               'error': str(error), 'error_type': type(error).__name__})
        raise
    return {'gate': 'one_chain_window_completed', 'production_ready': False,
            'completed_windows': index+1, 'receipt': str(directory/'charge-completed.json')}


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('action', choices=('inspect', 'run-next'))
    parser.add_argument('plan', type=Path)
    parser.add_argument('--max-steps', type=int, default=2000)
    parser.add_argument('--timeout', type=float, default=60.)
    args = parser.parse_args()
    try:
        report = inspect_plan(args.plan) if args.action == 'inspect' else run_next(args.plan, max_steps=args.max_steps, timeout=args.timeout)
    except (OSError, ValueError, KeyError, TypeError, OverflowError, subprocess.TimeoutExpired) as error:
        parser.exit(2, f'Charge chain failed: {error}\n')
    print(json.dumps(report, indent=2, sort_keys=True, allow_nan=False))


if __name__ == '__main__':
    main()
