from pathlib import Path
import time
import tempfile
import numpy as np
from qgpu.fep import stability, run_window, walk, pairs
from qgpu import Config, DistRestraint, Engine, SeqRestraint

from dataclasses import replace

FEP_DIR = Path(
    "~/code/qligfepv2-BenchmarkExperiments-shen/perturbations/cdk2/2.protein/FEP_1h1q_1oiu")


def pair_restraints(k):   # qligfep's heavy-atom pairs between the two ligands
    return tuple(DistRestraint(9033 + i, 9078 + i, 0.0, 0.1, k) for i in range(24))


def ligand_restraint(k):
    return (SeqRestraint(9033, 9127, k),)


def system(inputs, seed):  # everything that is the same in all 106 qligfep stages, plus eq1
    return Config(
        topology=inputs / "dualtop.top", fep=inputs / "FEP1.fep",
        seed=seed, initial_temperature=1,
        shell_radius=25, shell_force=10, lrf=True, lrf_cutoff=99, q_atom=99,
        solute_solute=10, solvent_solvent=10, solute_solvent=10,
        radial_force=60, polarisation=True, polarisation_force=20, non_bond=25,
        solute_constraint_algorithm="lincs", solvent_constraint_algorithm="lincs",
        separate_scaling=True, output=5, trajectory=0, lambdas=(0.5, 0.5),
        # eq1
        steps=10000, stepsize=0.1, temperature=1, bath_coupling=0.2,
        constrain_solvent=True, constrain_solute=False, constrain_hydrogens=False,
        seq_restraints=ligand_restraint(10.0), dist_restraints=pair_restraints(1.5),
    )


def create_fep_dir(replicate=1, temperature=298, out=None):
    """Returns the folder with the inputs and the output folder of this replicate."""
    fep_dir = FEP_DIR.expanduser()
    out = Path(out or fep_dir / "FEP1" /
               f"{temperature:g}" / str(replicate)).expanduser()
    out.mkdir(parents=True, exist_ok=True)
    return fep_dir, out


def create_engine(fep_dir, seed=2924, device=0, precision="spfp"):
    return Engine(system(fep_dir / "inputfiles", seed), device=device, precision=precision)


def run_equibrium(engine, out, temperature=298, scale=1.0):
    def steps(n):
        return max(1, round(n * scale))

    # eq1: fresh start, velocities from the seed
    engine.update_config(steps=steps(engine.config.steps))
    engine.run()

    # eq2..eq5: each continues from the one before; only what changes is passed
    # eq2
    engine.update_config(steps=steps(5000), stepsize=2.0, temperature=50, bath_coupling=2.0,
                         constrain_hydrogens=True)
    engine.run()

    # eq3
    engine.update_config(temperature=150, seq_restraints=ligand_restraint(5.0))
    engine.run()

    # eq4
    engine.update_config(temperature=275, seq_restraints=ligand_restraint(2.5))
    engine.run()

    # eq5
    part = out / "eq5.re.part"
    engine.update_config(steps=steps(50000), temperature=temperature, bath_coupling=10,
                         constrain_solute=True, seq_restraints=(), dist_restraints=pair_restraints(0.5),
                         final_file=str(part))
    engine.run()
    # otherwise the next run() overwrites it
    engine.update_config(final_file="")
    part.rename(out / "eq5.re")


def load_equilibrium(fep_dir, restart, temperature=298, device=0, precision="spfp"):
    """An engine that starts from eq5.re, with the settings eq5 ended with."""
    config = replace(
        system(fep_dir / "inputfiles", seed=None),
        # keep the velocities in the file
        restart=str(restart), initial_temperature=None,
        steps=5000, stepsize=2.0, temperature=temperature, bath_coupling=10,
        constrain_solute=True, constrain_hydrogens=True,
        seq_restraints=(), dist_restraints=pair_restraints(0.5),
    )
    return Engine(config, device=device, precision=precision)


def check_window(engine, out, lambdas=(0.5, 0.5), blocks=50, block=1000, energy=10,
                 limits=(0.5, 1.0, 1.5, 2.0)):
    """Run one window block by block and print what stability() says after each block."""
    rows, parts = [], []
    with tempfile.TemporaryDirectory(prefix="qgpu_check_") as scratch:
        for b in range(1, blocks + 1):
            engine.update_config(steps=block, lambdas=lambdas, energy=energy, output=block,
                                 energy_file=str(Path(scratch) / "block.en"))
            q = engine.run().q_energy
            parts.append(q[:, 0] - q[:, 1])
            g = np.concatenate(parts)
            rows.append(
                (b, parts[-1].mean(), [stability(g, limit) for limit in limits]))

    print(
        f"\nlambda {lambdas}: {block} steps per block, a frame every {energy} steps")
    print("block  block mean      mean   error   drift   stable at limit " +
          "  ".join(f"{x:g}" for x in limits))
    for b, block_mean, checks in rows:
        c = checks[0]
        marks = "    ".join("yes" if x.stable else " no" for x in checks)
        print(
            f"{b:5d}  {block_mean:10.2f}  {c.mean:8.2f}  {c.error:6.3f}  {c.drift:6.2f}   {marks}")

    print("\nblocks needed (at least 5 blocks, two stable checks in a row):")
    for i, limit in enumerate(limits):
        stop = next((b for b, _, checks in rows[1:]
                     if b >= 5 and checks[i].stable and rows[b - 2][2][i].stable), None)
        print(f"  limit {limit:g}: " + (
            f"{stop} blocks = {stop * block} steps" if stop else f"not within {blocks} blocks"))

    np.save(out / f"g_{round(lambdas[0] * 1000):04d}.npy", g)
    return g


def run(scale=1.0):
    fep_dir, out = create_fep_dir()
    for old in out.glob("md_*.en*"):
        old.unlink()

    restart = out / "eq5.re"
    if restart.exists():
        engine = load_equilibrium(fep_dir, restart)
    else:
        engine = create_engine(fep_dir)
        run_equibrium(engine, out, scale=scale)
    # check_window(engine, out)

    engine.update_config(output=1000)
    windows = walk(engine, out)
    for a, b, check in pairs(windows, engine.config.temperature):
        print(a.lambdas[0], b.lambdas[0], a.steps, round(
            check.forward, 2), round(check.backward, 2), check.fits)
    # window = run_window(engine, (0.5, 0.5), out / "md_0500_0500.en")

    # print(window.steps, window.stable, window.check)


if __name__ == "__main__":
    run()
