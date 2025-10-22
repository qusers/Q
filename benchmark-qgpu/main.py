import argparse
from benchmark_run import run 
    

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='Benchmark Qdyn runs (run/report).')

    parser.add_argument('--data_dir', type=str, help='Directory containing a single test case.')
    parser.add_argument('--bin', type=str, help='Path to the Qdyn GPU executable.')
    parser.add_argument('--max_processes', type=int, help='Max number of parallel processes to run.')

    args = parser.parse_args()

    run(args)