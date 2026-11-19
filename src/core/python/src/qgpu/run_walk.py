"""Run the adaptive lambda walk for one replicate, inside its run folder.

    python -m qgpu.run_walk [--template md_0500_0500.inp] [--dg-limit 0.07] ...

Everything about the system comes from the template .inp file: topology, FEP file,
temperature, cut-offs, restraints and the restart file to start from (eq5.re). The runner
has already filled in the temperature and the FEP file. Nothing is specific to one edge.

Written to the run folder:
    md_XXXX_YYYY.en   one energy file per window
    walk.tsv          one line per window: length, stability, overlap with the next window
    qfep.inp          rewritten to list the windows that were run
"""
import argparse
from dataclasses import replace
from pathlib import Path

from qgpu import Config, Engine
from qgpu.fep import _energy_file, pairs, walk


def read_config(template):
    """The settings of the template, without its own output files."""
    config = Config.from_inp(template)
    energy = config.energy
    # The walk writes the energy files itself; no trajectory and no restart files.
    # seed=None keeps the velocities of the restart file, as qdyn does for an md stage.
    config = replace(config, seed=None, initial_temperature=None,
                     output=1000, energy=0, energy_file="",
                     trajectory=0, trajectory_file="", final_file="")
    return config, energy


def write_qfep_inp(path, windows):
    """Rewrite qfep.inp for the windows that were run: the count and the file list."""
    path = Path(path)
    lines = path.read_text().splitlines()
    header = [line for line in lines[1:] if not line.startswith("md_")]
    # qfep.inp lists the files from lambda_1 = 1 down to 0
    files = [_energy_file(".", w.lambdas).name for w in reversed(windows)]
    path.write_text("\n".join([str(len(files))] + header + files) + "\n")


def write_table(path, windows, temperature):
    """walk.tsv: what the walk did, for later analysis."""
    overlaps = {a.lambdas[0]: check for a, _,
                check in pairs(windows, temperature)}
    with open(path, "w") as table:
        table.write(
            "lambda1\tlambda2\tsteps\tstable\tmean_gap\terror\tess_forward\tess_backward\n")
        for w in windows:
            # the last window has no next one
            check = overlaps.get(w.lambdas[0])
            ess = (f"{check.forward:.3f}",
                   f"{check.backward:.3f}") if check else ("", "")
            table.write(f"{w.lambdas[0]:.3f}\t{w.lambdas[1]:.3f}\t{w.steps}\t{int(w.stable)}\t"
                        f"{w.check.mean:.3f}\t{w.check.error:.3f}\t{ess[0]}\t{ess[1]}\n")


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    parser.add_argument("--template", default="md_0500_0500.inp",
                        help="a filled-in md .inp file of this replicate")
    parser.add_argument("--qfep-inp", default="qfep.inp")
    parser.add_argument("--ends", type=float, nargs=2, default=(0.999, 0.001),
                        metavar=("UP", "DOWN"), help="lambda_1 where the walk stops")
    parser.add_argument("--dg-limit", type=float, default=0.07,
                        help="largest error one window may add to dG, kcal/mol")
    parser.add_argument("--block", type=int, default=1000,
                        help="steps per block")
    parser.add_argument("--min-blocks", type=int, default=10)
    parser.add_argument("--max-blocks", type=int, default=50)
    parser.add_argument("--device", type=int, default=0)
    parser.add_argument("--precision", default="spfp")
    args = parser.parse_args(argv)

    config, energy = read_config(args.template)
    engine = Engine(config, device=args.device, precision=args.precision)
    windows = walk(engine, Path("."), ends=tuple(args.ends), dg_limit=args.dg_limit,
                   block=args.block, energy=energy,
                   min_blocks=args.min_blocks, max_blocks=args.max_blocks)

    write_table("walk.tsv", windows, config.temperature)
    write_qfep_inp(args.qfep_inp, windows)
    steps = sum(w.steps for w in windows)
    unstable = sum(not w.stable for w in windows)
    print(f"walk: {len(windows)} windows, {steps} steps, {unstable} not stable")


if __name__ == "__main__":
    main()
