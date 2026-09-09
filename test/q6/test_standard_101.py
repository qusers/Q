"""Standard array staging and endpoint-analysis contract."""
import importlib.util
from pathlib import Path
import shutil
import subprocess

import pytest

from QligFEP import charge_protocol as cp

ROOT = Path(__file__).resolve().parents[2]
spec = importlib.util.spec_from_file_location('standard_101', ROOT/'experiments/analytical-charge-validation/standard_101.py')
standard = importlib.util.module_from_spec(spec)
spec.loader.exec_module(standard)


@pytest.mark.parametrize('step,steps,wanted', [(2., 5000, 10000), (2., 50000, 100000), (.2, 5000, 5000)])
def test_adaptation_preserves_duration_restraints_and_restart(tmp_path, step, steps, wanted):
    src = tmp_path/'old.inp'
    src.write_text(f'''[MD]
steps {steps}
stepsize {step}
temperature T_VAR
shake_solvent on
shake_hydrogens on
[cut-offs]
q_atom 99
[sphere]
shell_radius 20
[solvent]
polarisation on
radial_force 60
polarisation_force 20
[intervals]
energy 10
output 25
trajectory 100
non_bond 25
[files]
topology dualtop.top
restart eq5.re
final md.re
energy md.en
trajectory md.dcd
[lambdas]
0.500 0.500
[distance_restraints]
1 2 0 .1 .5 0
[correction]
kernel 1
[trajectory_atoms]
not excluded
''')
    before = src.read_bytes()
    dest = tmp_path/'new.inp'
    dest.write_text(standard.adapt_input(src, 18.54))
    raw = cp.sections(dest)
    md = cp.keyed(raw['md'])
    assert int(md['steps']) == wanted
    assert float(md['stepsize'])*wanted == step*steps
    assert md['constraint_algorithm'] == 'shake shake' and md['lrf'] == 'off'
    assert set(cp.keyed(raw['cut-offs']).values()) == {'120'}
    assert cp.keyed(raw['files'])['restart'] == 'eq5.re'
    assert 'trajectory' not in cp.keyed(raw['files']) and 'correction' not in raw
    assert raw['distance_restraints'] == cp.sections(src)['distance_restraints']
    sol = cp.keyed(raw['solvent'])
    assert sol['polarization_adaptation'] == 'off' and sol['perstate_born_correction'] == 'on'
    assert sol['radius'] == '18.54' and src.read_bytes() == before


def test_qfep_excludes_exact_endpoint_files_only(tmp_path):
    src = tmp_path/'qfep.inp'
    names = ['md_1000_0000.en']+[f'md_{1000-i:04d}_{i:04d}.en' for i in range(1,100)]+['md_0000_1000.en']
    header = ['101', '2 0', '.592 100', '100', '100', '100', '0', '0', '1 0']
    src.write_text('\n'.join(header+names))
    trimmed = standard.trim_qfep(src).splitlines()
    assert trimmed[0] == '99' and trimmed[1:9] == header[1:]
    assert trimmed[9:] == names[1:-1]
    assert src.read_text().splitlines()[0] == '101'


def test_nonstandard_qfep_order_is_rejected(tmp_path):
    src = tmp_path/'qfep.inp'
    src.write_text('\n'.join(['101']+['0']*8+['md_0500_0500.en']*101))
    with pytest.raises(AssertionError):
        standard.trim_qfep(src)


def test_other_clusters_not_silently_selected(tmp_path):
    with pytest.raises(ValueError, match='SNELLIUS'):
        standard.stage(tmp_path, tmp_path/'new', '/project', 'LOCAL')


@pytest.mark.parametrize('target,edge', [
    ('cmet', 'FEP_CHEMBL3402742_23_CHEMBL3402744_300'),
    ('eg5', 'FEP_CHEMBL1085666_CHEMBL1089056')])
def test_private_target_fresh_start_when_available(tmp_path, target, edge):
    runtime = ROOT/'experiments/analytical-charge-validation/runtime'
    inputs = runtime/'standard-101-stage'/f'{target}-fwd/2.protein'/edge/'inputfiles'
    binary = runtime/'integration-20260908/build/src/q6/qdyn'
    if not inputs.exists() or not binary.exists():
        pytest.skip('Optional private prepared targets and retained local build are absent')
    for name in ('dualtop.top', 'FEP1.fep'):
        shutil.copyfile(inputs/name, tmp_path/name)
    raw = cp.sections(inputs/'eq1.inp')
    md = cp.keyed(raw['md'])
    md.update(steps='100', random_seed='2924')
    raw['md'] = [line.split() for line in (k+' '+v for k, v in md.items())]
    raw['files'] = [[k, 'FEP1.fep' if k == 'fep' else v] for k, v in cp.keyed(raw['files']).items()]
    inp = tmp_path/'eq1.inp'
    inp.write_text(standard.serialize(raw))
    log = tmp_path/'eq1.log'
    with log.open('w') as stream:
        subprocess.run([str(binary), inp.name], cwd=tmp_path, stdout=stream, stderr=subprocess.STDOUT,
                       timeout=120, check=True)
    standard.check(inp, log)
