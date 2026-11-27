from pathlib import Path
import time

from qgpu import Config, Engine

INPUT = Path("~/code/qligfepv2-BenchmarkExperiments-shen/perturbations/experiment/cdk2_test/eq5.inp")

config = Config.from_inp(INPUT.expanduser())

start = time.perf_counter()
engine = Engine(config, device=0)
engine.run()

print(f"Completed {config.steps} steps of {engine.n_atoms} atoms in {time.perf_counter() - start:.2f} seconds")

