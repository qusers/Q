import argparse


def positive_int(value):
    parsed = int(value)
    if parsed < 1:
        raise argparse.ArgumentTypeError("must be >= 1")
    return parsed
    

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Benchmark native .inp Qdyn runs.")

    parser.add_argument("--input", required=True, help="Path to the Qdyn .inp file to benchmark.")
    parser.add_argument("--bin", required=True, help="Path to the Qdyn executable.")
    execution_mode = parser.add_mutually_exclusive_group(required=True)
    execution_mode.add_argument(
        "--max-processes",
        "--max_processes",
        dest="max_processes",
        type=positive_int,
        help="Sweep from 1 through this many independent Qdyn processes.",
    )
    execution_mode.add_argument(
        '--concurrency',
        type=positive_int,
        nargs='+',
        help='Specific parallel process counts to run, e.g. --concurrency 1 2 3 4 5 10 15 20.',
    )

    args = parser.parse_args()

    from benchmark_run import run

    run(args)
