"""Run replica exchange over fixed lambda windows for one replicate (prototype).

    python -m qgpu.run_hrex --cycles 500 [--windows 16] [--workers 16] [--rest-min 0.5] [--no-exchange]
    python -m qgpu.run_hrex --lambdas 0 0.05 0.1 ... 1 --cycles 500

Run inside a replicate folder after eq5, like qgpu.run_walk. The windows are given, either
as a number of evenly spaced ones (0 and 1 included; 16 by default) or as their lambda_1 values. Each gets
a starting configuration from a short run (see hrex.prepare); replica exchange then runs
`cycles` blocks in every window.

Written to the run folder:
    md_XXXX_YYYY.en        one energy file per window
    md_XXXX_YYYY.en.fhe    foreign work of each window to itself and its neighbours
    qfep.inp               rewritten to list the windows (qfep is only right without soft-core)
    hrex.tsv               per window: lambda, REST2 scale, acceptance with the next window, and the free
                           energy step to it (BAR on the foreign work, kcal/mol)
    hrex_dg.txt            free energy of state 1 -> state 2: all frames, and each half
    hrex.npz               replica (cycles+1, n), coords (cycles, n, n_saved, 3), lambdas,
                           round_trips, save_atoms, dg_steps, work_up and work_down
                           (see hrex.neighbour_work)
"""
import argparse
import os
from pathlib import Path

import numpy as np

from qgpu import Engine
from qgpu.hrex import Workers, even_windows, rest_ladder, window_pairs, free_energy, hrex, neighbour_work, prepare, round_trips
from qgpu.run_walk import read_config, write_qfep_inp


def atom_range(text):
    """'first:last[,first:last...]', 1-based and inclusive as in the topology -> 0-based indices."""
    parts = []
    for item in text.split(","):
        first, last = (int(x) for x in item.split(":"))
        if first < 1 or last < first:
            raise argparse.ArgumentTypeError("use first:last with 1 <= first <= last")
        parts.append(np.arange(first - 1, last))
    return np.concatenate(parts)


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
    which = parser.add_mutually_exclusive_group()
    which.add_argument("--windows", type=int, default=16,
                       help="number of evenly spaced windows, 0 and 1 included")
    which.add_argument("--lambdas", type=float, nargs="+", metavar="L1", help="lambda_1 of every window")
    parser.add_argument("--rest-min", type=float, default=1.0,
                        help="REST2 scale of the Q atoms in the middle window (T0 / highest effective "
                             "temperature); the ends stay at 1. 1 = no REST (see hrex.rest_ladder)")
    parser.add_argument("--prepare-steps", type=int, default=5000,
                        help="steps each window runs before replica exchange starts")
    parser.add_argument("--cycles", type=int, required=True, help="blocks per window")
    parser.add_argument("--block", type=int, default=500, help="steps between swap attempts")
    parser.add_argument("--no-exchange", action="store_true", help="control: same schedule, no swaps")
    parser.add_argument("--save-atoms", type=atom_range, default=None, metavar="FIRST:LAST[,...]",
                        help="keep the coordinates of these atoms after every block")
    parser.add_argument("--save-every", type=int, default=1,
                        help="keep the coordinates of --save-atoms after every N-th block only")
    parser.add_argument("--work-dir", default=os.environ.get("TMPDIR"),
                        help="folder for the files that grow during the run (default: $TMPDIR); "
                             "they move to the run folder at the end")
    parser.add_argument("--workers", type=int, default=1,
                        help="processes that run the windows of a cycle at the same time, on one GPU")
    parser.add_argument("--seed", type=int, default=0)
    parser.add_argument("--device", type=int, default=0)
    parser.add_argument("--precision", default="spfp")
    args = parser.parse_args(argv)

    config, energy = read_config(args.template)
    lambda1 = window_pairs(args.lambdas if args.lambdas else even_windows(args.windows))

    def run(engine, prepare_windows):
        windows = prepare_windows(lambda1, steps=args.prepare_steps, scales=rest_ladder(len(lambda1), args.rest_min))
        print(f"prepared {len(windows)} windows, {args.prepare_steps} steps each", flush=True)
        return hrex(engine, windows, Path("."), cycles=args.cycles, block=args.block, energy=energy,
                    exchange=not args.no_exchange, save_atoms=args.save_atoms, seed=args.seed,
                    work_dir=args.work_dir, save_every=args.save_every)

    if args.workers > 1:
        with Workers(config, min(args.workers, len(lambda1)), device=args.device, precision=args.precision) as workers:
            result = run(workers, workers.prepare)
    else:
        engine = Engine(config, device=args.device, precision=args.precision)
        result = run(engine, lambda l, steps, scales: prepare(engine, l, steps=steps, scales=scales))

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
        table.write("lambda1\tlambda2\trest_scale\ttried\taccepted\tacceptance_with_next\tdF_to_next\n")
        for k, l in enumerate(result.lambdas):
            pair = (f"{result.tried[k]}\t{result.accepted[k]}\t{result.acceptance[k]:.3f}\t{dg_steps[k]:.4f}"
                    if k < len(result.tried) else "\t\t\t")
            table.write(f"{l[0]:.3f}\t{l[1]:.3f}\t{result.scales[k]:.4f}\t{pair}\n")
    Path("hrex_dg.txt").write_text(
        f"dG_bar_foreign {-dg_steps.sum():.4f}\nfirst_half {first:.4f}\nsecond_half {second:.4f}\n")
    print(f"dG (state 1 -> 2, BAR on foreign work): {-dg_steps.sum():.3f} kcal/mol "
          f"(first half {first:.3f}, second half {second:.3f})")
    work_up, work_down = neighbour_work(result)
    np.savez_compressed("hrex.npz", replica=result.replica, coords=result.coords,
                        lambdas=np.array(result.lambdas), round_trips=trips, dg_steps=dg_steps,
                        scales=np.array(result.scales), save_every=args.save_every,
                        work_up=work_up.astype(np.float32), work_down=work_down.astype(np.float32),
                        save_atoms=np.array([] if args.save_atoms is None else args.save_atoms))
    steps = len(result.lambdas) * args.cycles * args.block
    print(f"hrex: {len(result.lambdas)} windows, {steps} steps, "
          f"mean acceptance {result.acceptance.mean():.2f}, round trips {int(trips.sum())}")


if __name__ == "__main__":
    main()
