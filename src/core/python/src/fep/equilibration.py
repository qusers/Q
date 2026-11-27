"""Heating and relaxation of one end state before the windows are made."""
import os
import time
from dataclasses import dataclass, replace
from pathlib import Path

from fep.log import log


@dataclass(frozen=True)
class Stage:
    steps: int
    stepsize: float          # fs
    temperature: float       # K
    bath_coupling: float
    constrain_hydrogens: bool
    constrain_solute: bool
    k_ligand: float          # restraint on the ligand atoms
    to_center: int           # 1: restrain the ligand's centre instead of each atom
    k_pair: float            # restraints between the two ligand copies


def heating(temperature):
    """
    Heating and relaxation in five stages, as eq1-eq5 of QligFEP
    Updates: Always set constrain_solute false
    """
    return (
        Stage(10000, 0.1, 1, 0.2, False, False, 10.0, 0, 1.5),
        Stage(5000, 2.0, 50, 2.0, True, False, 10.0, 0, 1.5),
        Stage(5000, 2.0, 150, 2.0, True, False, 5.0, 0, 1.5),
        Stage(5000, 2.0, 275, 2.0, True, False, 2.5, 0, 1.5),
        Stage(50000, 2.0, temperature, 10.0, True, False, 1.0, 1, 0.5),
    )


def equilibrate(engine, config, lambdas, restart):
    """Heat and relax the system at `lambdas`; returns the final state.

    The final state is also written to the Q restart file `restart`. If that file is already
    there, it is loaded and nothing is run. Otherwise the run starts from the engine's current
    state. `config` is the template: its restraints give the atoms, the stages give the strengths.
    """
    restart = Path(restart)
    if restart.exists():
        log.info("end state lambdas %s: loaded from %s, no heating", lambdas, restart)
        engine.load_restart(restart)
        return engine.get_state()
    # The engine also writes the file during the run: use another name until the run is over,
    # so that a killed job leaves no file that looks finished.
    unfinished = restart.with_name(restart.name + ".unfinished")
    stages = heating(config.temperature)
    started = time.time()
    for number, stage in enumerate(stages, 1):
        log.info("end state lambdas %s: heating stage %d of %d, %d steps of %g fs at %g K",
                 lambdas, number, len(stages), stage.steps, stage.stepsize, stage.temperature)
        engine.update_config(
            steps=stage.steps, stepsize=stage.stepsize, temperature=stage.temperature,
            bath_coupling=stage.bath_coupling, constrain_hydrogens=stage.constrain_hydrogens,
            constrain_solute=stage.constrain_solute, lambdas=lambdas,
            seq_restraints=[replace(r, k=stage.k_ligand, to_center=stage.to_center)
                            for r in config.seq_restraints],
            dist_restraints=[replace(r, k=stage.k_pair) for r in config.dist_restraints],
            final_file=str(unfinished) if stage is stages[-1] else "")
        result = engine.run()
        log.info("    done: potential energy %.1f kcal/mol", result.final_upot)
    engine.update_config(final_file="")
    os.replace(unfinished, restart)
    log.info("end state lambdas %s: heated in %.0f s, written to %s", lambdas, time.time() - started, restart)
    return engine.get_state()
