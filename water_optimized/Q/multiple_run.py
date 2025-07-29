import subprocess
import time
import multiprocessing
import os
import argparse

def run_program(program_index, program_command, output_file):
    """
    Run a program and write its output to a specified file, recording the execution time.
    """
    # Expand the ~ symbol in paths
    program_command = os.path.expanduser(program_command)
    output_file = os.path.expanduser(output_file)

    # Record the start time for each process
    start_time = time.time()

    try:
        process = subprocess.Popen(
            program_command, stdout=open(output_file, "w"), stderr=subprocess.PIPE, shell=True
        )
        # Wait for the process to complete
        process.wait()

        end_time = time.time()
        execution_time = end_time - start_time
        print(f"Process {program_index + 1} (output to {output_file}) finished, execution time: {execution_time:.2f} seconds")

    except Exception as e:
        print(f"Process {program_index + 1} failed with error: {e}")


def main(program_command, process_num):
    # Define the output files for each process
    outputs = [f"output_{i+1}.log" for i in range(process_num)]

    # Open a pool of processes
    with multiprocessing.Pool(processes=process_num) as pool:
        results = [
            pool.apply_async(run_program, (i, program_command, outputs[i]))
            for i in range(process_num)
        ]

        # 等待所有进程完成
        for result in results:
            result.get()

    print("All processes finished.")


if __name__ == "__main__":
    # Use argparse to handle command line arguments
    parser = argparse.ArgumentParser(description='Run multiple instances of a program.')
    parser.add_argument('process_num', type=int, help='Number of processes to run simultaneously')
    parser.add_argument('program_command', type=str, help='Command to run the program (with path)')

    args = parser.parse_args()

    main(args.program_command, args.process_num)
