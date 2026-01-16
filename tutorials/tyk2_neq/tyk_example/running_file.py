import os
import sys
import subprocess
import time
if len(sys.argv) != 5:
    sys.stderr.write("Usage: python running_file.py [folder] [input] [output] [qdyn_directory]\n")
    sys.exit(1)

folder, input_name, output_name, qdyn_dir = sys.argv[1:5]

try:
    print(f'cd {folder}')
    os.chdir(folder)
    time.sleep(0.5)
    print(os.popen('pwd').read())
    cmd = f"{qdyn_dir}/qdyn_neq {input_name}.inp > {output_name}.log &"
    os.popen(cmd).read()
except Exception as e:
    sys.stderr.write(f"Failed to change directory to {folder}: {e}\n")
    sys.exit(2)





