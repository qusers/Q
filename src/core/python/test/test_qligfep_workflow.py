"""The qligfep FEP workflow (FEP_1h1q_1oiu): every replicate in its own process.

Each replicate is run_one_replicate() from test_qligfep_one_replicate.py.
Writes md_XXXX_YYYY.en files for qfep.
"""

from concurrent.futures import ProcessPoolExecutor
import multiprocessing as mp
from pathlib import Path

from test_qligfep_one_replicate import run_one_replicate

SEEDS = [2924, 25360, 21448, 14380, 14188, 28133, 2816, 22850, 6601, 3085]


def run_workflow(replicates=len(SEEDS), gpus=1, temperature=298, precision="spfp", scale=1.0, out=None):
    """Replicates 1..`replicates`, one process each, spread over `gpus` GPUs.

    `out`: replicate i goes to <out>/<i> (default <fep-dir>/FEP1/<T>/<i>); `scale` shortens every stage.
    """
    with ProcessPoolExecutor(replicates, mp_context=mp.get_context("spawn")) as pool:
        jobs = [pool.submit(run_one_replicate, i, SEEDS[i - 1], temperature, precision, (i - 1) % gpus, scale,
                            Path(out) / str(i) if out else None)
                for i in range(1, replicates + 1)]
        for job in jobs:
            job.result()  # re-raises an error from that replicate


if __name__ == "__main__":
    run_workflow()
