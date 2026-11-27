"""Engine.update_config(lambdas=...) must give the same dynamics as an engine built at that lambda.

Needs the GPU and eq5.re from test_fep.py. Both engines start from the same restart file and run
60 steps at (0.9, 0.1); over so few steps two identical engines agree exactly.
"""

from dataclasses import replace

import numpy as np

from qgpu import Engine
from test_fep import create_fep_dir, load_equilibrium

LAMBDAS = (0.9, 0.1)


def short_run(engine, energy_file):
    engine.update_config(lambdas=LAMBDAS, steps=60, energy=10,
                         output=60, energy_file=str(energy_file))
    return engine.run().q_energy


def test_update_config_lambdas_matches_an_engine_built_at_that_lambda(tmp_path):
    fep_dir, out = create_fep_dir()

    # built at (0.5, 0.5)
    engine = load_equilibrium(fep_dir, out / "eq5.re")
    built_config = replace(engine.config, lambdas=LAMBDAS)
    updated = short_run(engine, tmp_path / "updated.en")
    del engine

    built = short_run(Engine(built_config, device=0,
                      precision="spfp"), tmp_path / "built.en")

    assert np.abs(
        updated - built).max() < 1e-3, f"updated\n{updated}\nbuilt\n{built}"
