"""One replicate of the qligfep FEP workflow (FEP_1h1q_1oiu), in this process, without .inp files.

eq1..eq5, then the lambda windows from 0.5 toward 1.0 and toward 0.0. Writes md_XXXX_YYYY.en
files for qfep. test_qligfep_workflow.py runs this for every replicate in parallel.
"""

from pathlib import Path
import time

from qgpu import Config, DistRestraint, Engine, SeqRestraint

FEP_DIR = Path(
    "~/code/qligfepv2-BenchmarkExperiments-shen/perturbations/cdk2/2.protein/FEP_1h1q_1oiu")


def sigmoidal_lambdas(windows=100):
    """qligfep's sigmoidal schedule (IO.get_lambdas) as (lambda1, lambda2) pairs, sorted by lambda1.

    Both values are rounded separately, as in qligfep, so a pair can sum to 1.001 (0.988, 0.013).
    """
    k, half = -1.1, windows // 2
    sigmoid = [k * (i / half) / (k - i / half + 1.0) for i in range(half + 1)]
    upper = [round(0.5 * (s + 1), 3) for s in sigmoid]   # 0.5 -> 1.0
    lower = [round(0.5 * (1 - s), 3) for s in sigmoid]   # 0.5 -> 0.0
    return sorted(set(zip(upper, lower)) | set(zip(lower, upper)))


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


def run_one_replicate(replicate=1, seed=2924, temperature=298, precision="spfp", device=0, scale=1.0, out=None):
    """`seed` (> 0) sets the eq1 velocities; `replicate` names the output folder; `scale` shortens every stage (for testing)."""
    def steps(n):
        return max(1, round(n * scale))

    fep_dir = FEP_DIR.expanduser()
    out = Path(out or fep_dir / "FEP1" /
               f"{temperature:g}" / str(replicate)).expanduser()
    out.mkdir(parents=True, exist_ok=True)
    t0 = time.perf_counter()

    # eq1: fresh start, velocities from the seed
    config = system(fep_dir / "inputfiles", seed)
    engine = Engine(config, device=device, precision=precision)
    engine.update_config(steps=steps(config.steps))
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
    engine.update_config(steps=steps(50000), temperature=temperature, bath_coupling=10,
                         constrain_solute=True, seq_restraints=(), dist_restraints=pair_restraints(0.5))
    engine.run()

    # FEP windows: one run per lambda pair, each writing md_XXXX_YYYY.en
    def window(l1, l2):
        name = f"md_{round(l1 * 1000):04d}_{round(l2 * 1000):04d}.en"
        engine.update_config(steps=steps(5000), lambdas=(l1, l2), energy=10, energy_file=str(out / name))
        engine.run()

    lambdas = sigmoidal_lambdas()

    # 1. the midpoint 0.5/0.5; both branches start from its final state
    window(0.5, 0.5)
    midpoint = engine.get_state()

    # 2. toward state 1: 0.592 -> 1.0, continuing from the midpoint
    for l1, l2 in [p for p in lambdas if p[0] > 0.5]:
        window(l1, l2)

    # 3. toward state 2: 0.408 -> 0.0, starting again from the midpoint
    engine.set_state(midpoint)
    for l1, l2 in [p for p in reversed(lambdas) if p[0] < 0.5]:
        window(l1, l2)

    print(
        f"replicate {replicate}: {len(list(out.glob('md_*.en')))} windows in {time.perf_counter() - t0:.0f} s -> {out}")


if __name__ == "__main__":
    run_one_replicate()
