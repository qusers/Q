import argparse


def positive_int(value):
    parsed = int(value)
    if parsed < 1:
        raise argparse.ArgumentTypeError("must be >= 1")
    return parsed
    

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='Benchmark Qdyn runs (run/report).')

    parser.add_argument('--data_dir', type=str, help='Directory containing a single test case.')
    parser.add_argument('--bin', type=str, help='Path to the Qdyn GPU executable.')
    parser.add_argument('--max_processes', type=positive_int, help='Max number of parallel processes to run.')
    parser.add_argument(
        '--concurrency',
        type=positive_int,
        nargs='+',
        help='Specific parallel process counts to run, e.g. --concurrency 1 2 3 4 5 10 15 20.',
    )

    args = parser.parse_args()

    from benchmark_run import run

    run(args)
