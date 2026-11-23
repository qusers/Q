"""Scan the lambda windows of one replicate, then run replica exchange over them (prototype).

    python -m qgpu.run_hrex --cycles 500 [--block 500] [--no-exchange] [--save-atoms 9033:9131]

Run inside a replicate folder after eq5, like qgpu.run_walk. The scan is the adaptive
walk with its usual settings; it places the windows and gives each one a starting
configuration. Replica exchange then runs `cycles` blocks in every window.

Written to the run folder:
    scan/md_XXXX_YYYY.en   energies of the scan (not used for the result)
    md_XXXX_YYYY.en        one energy file per window, replica exchange only
    qfep.inp               rewritten to list the windows
    md_XXXX_YYYY.en.fhe    foreign work of each window to itself and its neighbours
    hrex.tsv               per window: lambda, acceptance with the next window, and the free
                           energy step to it (BAR on the foreign work, kcal/mol)
    hrex_dg.txt            free energy of state 1 -> state 2: all frames, and each half
    hrex.npz               replica (cycles+1, n), coords (cycles, n, n_saved, 3), lambdas,
                           round_trips, save_atoms, dg_steps, work_up and work_down
                           (see hrex.neighbour_work)
"""
import argparse
from pathlib import Path

import numpy as np

from qgpu import Engine
from qgpu.fep import walk
from qgpu.hrex import free_energy, hrex, neighbour_work, round_trips
from qgpu.run_walk import read_config, write_qfep_inp


def atom_range(text):
    """'first:last', 1-based and inclusive as in the topology -> 0-based indices."""
    first, last = (int(x) for x in text.split(":"))
    if first < 1 or last < first:
        raise argparse.ArgumentTypeError("use first:last with 1 <= first <= last")
    return np.arange(first - 1, last)


def _first(result, frames):
    """The result with only the first `frames` frames of every window."""
    from dataclasses import replace
    return replace(result, foreign_work=[w[:frames] for w in result.foreign_work])


class _Placed:
    def __init__(self, lambdas):
        self.lambdas = lambdas


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    parser.add_argument("--template", default="md_0500_0500.inp")
    parser.add_argument("--qfep-inp", default="qfep.inp")
    parser.add_argument("--ends", type=float, nargs=2, default=(0.999, 0.001), metavar=("UP", "DOWN"))
    parser.add_argument("--dg-limit", type=float, default=0.07, help="scan: see qgpu.run_walk")
    parser.add_argument("--scan-block", type=int, default=1000)
    parser.add_argument("--scan-min-blocks", type=int, default=10)
    parser.add_argument("--scan-max-blocks", type=int, default=50)
    parser.add_argument("--cycles", type=int, required=True, help="blocks per window")
    parser.add_argument("--block", type=int, default=500, help="steps between swap attempts")
    parser.add_argument("--no-exchange", action="store_true", help="control: same schedule, no swaps")
    parser.add_argument("--save-atoms", type=atom_range, default=None, metavar="FIRST:LAST",
                        help="keep the coordinates of these atoms after every block")
    parser.add_argument("--seed", type=int, default=0)
    parser.add_argument("--device", type=int, default=0)
    parser.add_argument("--precision", default="spfp")
    args = parser.parse_args(argv)

    config, energy = read_config(args.template)
    engine = Engine(config, device=args.device, precision=args.precision)

    scan_dir = Path("scan")
    scan_dir.mkdir(exist_ok=True)
    windows = walk(engine, scan_dir, ends=tuple(args.ends), dg_limit=args.dg_limit,
                   block=args.scan_block, energy=energy,
                   min_blocks=args.scan_min_blocks, max_blocks=args.scan_max_blocks)
    print(f"scan: {len(windows)} windows, {sum(w.steps for w in windows)} steps", flush=True)

    result = hrex(engine, windows, Path("."), cycles=args.cycles, block=args.block, energy=energy,
                  exchange=not args.no_exchange, save_atoms=args.save_atoms, seed=args.seed)

    write_qfep_inp(args.qfep_inp, [_Placed(l) for l in result.lambdas])
    trips = round_trips(result.replica)
    # BAR on the foreign work: right with and without soft-core (qfep is only right without)
    temperature = config.temperature
    dg_steps = free_energy(result, temperature)
    frames = min(len(w) for w in result.foreign_work)
    half = frames // 2
    first = -free_energy(_first(result, half), temperature).sum()
    second = -free_energy(result, temperature, skip=half).sum()
    with open("hrex.tsv", "w") as table:
        table.write("lambda1\tlambda2\ttried\taccepted\tacceptance_with_next\tdF_to_next\n")
        for k, l in enumerate(result.lambdas):
            pair = (f"{result.tried[k]}\t{result.accepted[k]}\t{result.acceptance[k]:.3f}\t{dg_steps[k]:.4f}"
                    if k < len(result.tried) else "\t\t\t")
            table.write(f"{l[0]:.3f}\t{l[1]:.3f}\t{pair}\n")
    Path("hrex_dg.txt").write_text(
        f"dG_bar_foreign {-dg_steps.sum():.4f}\nfirst_half {first:.4f}\nsecond_half {second:.4f}\n")
    print(f"dG (state 1 -> 2, BAR on foreign work): {-dg_steps.sum():.3f} kcal/mol "
          f"(first half {first:.3f}, second half {second:.3f})")
    work_up, work_down = neighbour_work(result)
    np.savez_compressed("hrex.npz", replica=result.replica, coords=result.coords,
                        lambdas=np.array(result.lambdas), round_trips=trips, dg_steps=dg_steps,
                        work_up=work_up.astype(np.float32), work_down=work_down.astype(np.float32),
                        save_atoms=np.array([] if args.save_atoms is None else args.save_atoms))
    steps = len(result.lambdas) * args.cycles * args.block
    print(f"hrex: {len(result.lambdas)} windows, {steps} steps, "
          f"mean acceptance {result.acceptance.mean():.2f}, round trips {int(trips.sum())}")


if __name__ == "__main__":
    main()
