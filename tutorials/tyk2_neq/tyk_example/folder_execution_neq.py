import random
import os
import re
import time
import sys

#this is a script that executes NEQ simulations
#it accepts the name of the folder to run the simulations, and the sim type (protein or water)
#the path to the Q executable and the stepsize needs to be written here
# also accepts failure conditions like timeouts in seconds

path_to_q = "/gpfs/home5/vprypoten/Qengine_neq/q6/bin/q6/qdyn_neq"
stepsize = 2.0  # Time step for the simulation

new_seed = True  # Set to False to keep the same seed

#different wait times for different stages. If it exceeds this, the script dies
eq_wait_timer = 600
eq6_wait_timer = 600
eq6_generic_wait_timer = 300
neq_wait_timer = 300






#default parameters if the parameter file is missing
number_of_reps = 5          #number of NEQ simulations per independent realization
number_of_steps_eq5 = 5000  # Number of steps for eq6_0 and eq6_1
number_of_steps_eq6 = 50  # Number of steps for eq6_0 and eq6_1
number_of_steps_neq = 50000  # Number of steps for neq_0 and neq_1
L_parameter = 8.0  # Default L parameter for sigmoidal scaling
lines_to_add = ['[lambda_scaling] ','scaling_parameter          sigmoidal', f'L_sigmoid        {L_parameter}']

parameter_file = "parameters.txt"  # File containing parameters for each realization
if not os.path.exists(parameter_file):
    print(f"Parameter file {parameter_file} not found. Using default parameters.")
else:
    with open(parameter_file, 'r') as f:
        lines = f.readlines()
    
    for line in lines:
        match_reps = re.match(r'number_of_reps\s*=\s*(\d+)', line)
        match_steps_eq5 = re.match(r'number_of_steps_eq5\s*=\s*(\d+)', line)
        match_steps_eq6 = re.match(r'number_of_steps_eq6\s*=\s*(\d+)', line)
        match_steps_neq = re.match(r'number_of_steps_neq\s*=\s*(\d+)', line)
        match_L_parameter = re.match(r'L_parameter\s*=\s*([0-9.]+)', line)
        elif match_reps:
            number_of_reps = int(match_reps.group(1))
        elif match_steps_eq5:
            number_of_steps_eq5 = int(match_steps_eq5.group(1))
        elif match_steps_eq6:
            number_of_steps_eq6 = int(match_steps_eq6.group(1))
        elif match_steps_neq:
            number_of_steps_neq = int(match_steps_neq.group(1))
        elif match_L_parameter:
            L_parameter = float(match_L_parameter.group(1))
            lines_to_add = ['[lambda_scaling] ','scaling_parameter          sigmoidal', f'L_sigmoid        {L_parameter}']
    print("Parameters loaded from file.")
    

# Example usage:
# replacements = {"key1=": "new_value1", "key2=": "new_value2"}
# copy_and_modify_file("water/subfolder/file.txt", "neq_water/subfolder/file.txt", replacements)
def copy_and_modify_file(input_path, output_path, replacements):
    """
    Reads input_path, replaces lines starting with keys in replacements,
    and writes the result to output_path.
    """
    with open(input_path, 'r') as infile, open(output_path, 'w') as outfile:
        for line in infile:
            replaced = False
            for key, new_value in replacements.items():
                if line.startswith(key):
                    outfile.write(f"{key}{new_value}\n")
                    replaced = True
                    break
            if not replaced:
                outfile.write(line)
                
def modify_file(input_path, inline, outline):
    """
    Reads input_path, replaces lines starting with inline,
    and writes the result back to input_path.
    """
    with open(input_path, 'r') as infile:
        lines = infile.readlines()
    with open(input_path, 'w') as outfile:
        for line in lines:
            if line.startswith(inline):
                outfile.write(outline + "\n")
            else:
                outfile.write(line)


def append_lines_to_file(file_path, lines_to_add):
    """
    Appends each line in lines_to_add to the end of file_path.
    lines_to_add should be a list of strings.
    """
    with open(file_path, 'a') as f:
        for line in lines_to_add:
            f.write(f"{line}\n")

def log_file_reader(log_path):
    """
    reads a log file, and checks if it is completed, returns:
    0 if incomplete, 1 if normal termination, 2 if error termination
    """
    if not os.path.exists(log_path):
        return 0
    with open(log_path, 'r') as f:
        log_content = f.read()
    f.close()
    if "terminated normally" in log_content:
        return 1
    elif "!!!!!!!!!!!!!!!!" in log_content:
        return 2
    else:
        return 0
#python running_file.py [folder] [input] [output] [qdyn_directory]

def script_prep_function(subfolder,sim_type, rep):


    if new_seed:
        random_seed = random.randint(1, 1000000)
    else:
        random_seed = 123456  # Fixed seed for reproducibility
    print(f"Preparing realization {rep} with seed {random_seed}")
    output_subfolder = './neq_'+sim_type +'/'+ subfolder+'_rep'+str(rep)
    # Create the output folder if it doesn't exist
    os.makedirs(output_subfolder, exist_ok=True)
    
    input_path = os.path.join(sim_type, subfolder, "inputfiles")
    #copy the topology (dualtop.top), FEP (FEP1.fep) to the output_subfolder
    os.system(f'cp {os.path.join(input_path, "dualtop.top")} {output_subfolder}/dualtop.top')
    os.system(f'cp {os.path.join(input_path, "FEP1.fep")} {output_subfolder}/FEP1.fep')
    time.sleep(0.2)  # Ensure the copy operation completes
    

    #replace EQ1 through EQ4 files
    for i in range(1, 5):
        input_file = os.path.join(input_path, f"eq{i}.inp")
        output_file = os.path.join(output_subfolder, f"eq{i}.inp")
        replacements = {
            "random_seed       ": f" {random_seed}",
            "fep               ": " FEP1.fep",
        }
        copy_and_modify_file(input_file, output_file, replacements)
        print(f"Processed: {output_file}")
        time.sleep(0.2)  # Ensure the copy operation completes
    #now do eq5, and set temperature
    input_file = os.path.join(input_path, "eq5.inp")
    output_file = os.path.join(output_subfolder, "eq5.inp")
    replacements = {

        "fep              ": " FEP1.fep",
        "temperature      ": " 300",
        "steps             ": f" {number_of_steps_eq5}",

    }
    copy_and_modify_file(input_file, output_file, replacements)
    time.sleep(0.2)  # Ensure the copy operation completes
    print(f"Processed: {output_file}")
    #now prepare for the eq6_0_0 through eq_6_0_n and same for eq6_1
    #first we make a default eq6_0 and eq6_1 file
    input_file = os.path.join(input_path, "eq5.inp")
    output_file = os.path.join(output_subfolder, "eq6_0_default.inp")
    replacements = {
        "fep               ": " FEP1.fep",
        "stepsize          ": f" {stepsize}",
        "temperature       ": " 300",
        "output            ": " 10",
        "steps             ": f" {number_of_steps_eq5}"
    }
    copy_and_modify_file(input_file, output_file, replacements)
    time.sleep(0.2)  # Ensure the copy operation completes
    modify_file(output_file, '0.500', '0.0   1.0')
    time.sleep(0.2)  # Ensure the copy operation completes
    print(f"Processed: {output_file}")

    output_file = os.path.join(output_subfolder, "eq6_1_default.inp")
    replacements = {
        "fep               ": " FEP1.fep",
        "stepsize          ": f" {stepsize}",
        "temperature       ": " 300",
        "output            ": " 10",
        "steps             ": f" {number_of_steps_eq5}"

    }       
    copy_and_modify_file(input_file, output_file, replacements)
    time.sleep(0.2)  # Ensure the copy operation completes
    modify_file(output_file, '0.500', '1.0   0.0')
    time.sleep(0.2)  # Ensure the copy operation completes
    print(f"Processed: {output_file}")
        

def clean_up_function(subfolder, sim_type, index_of_run, rep):
    #this function cleans up space, by removing .re, .log files from eq6_?_n simulations that are no longer needed
    #it will also remove slurm*out files
    cleaning_directory = './neq_'+sim_type +'/'+ subfolder+'_rep'
    
    output_subfolder = cleaning_directory+str(rep)
    print(f"Cleaning up: {output_subfolder}")
    clean_from = index_of_run - 20  #clean files that are 20 below the index
    clean_to = index_of_run - 10  #clean files that are 10 below the index

    #remove neq.inp, eq6_0_n and eq6_1_n files that are no longer needed. they are 10 below the index
    for i in range(2):  # 0 and 1
        for cleaning_variable in range(clean_from,clean_to):
                    cleaning_name_eq = f'eq6_[0-1]_rep{cleaning_variable}.*'
                    cleaning_name_eqlog = f'eq6_[0-1]_{cleaning_variable}.log'
                    cleaning_name_neq = f'neq_[0-1]_rep{cleaning_variable}.inp'
                    os.popen(f'rm {output_subfolder}/{cleaning_name_eq}').read()
                    os.popen(f'rm {output_subfolder}/{cleaning_name_neq}').read()
                    os.popen(f'rm {output_subfolder}/{cleaning_name_eqlog}').read()
    os.popen(f'rm {output_subfolder}/slurm*out').read()
    print('cleaning done')

#now, we execute the simulations
#we will use functions so we can call on them only if the simulations were not previously completed
def equilibration_function(subfolder, sim_type, rep):
    

    output_subfolder = './neq_'+sim_type +'/'+ subfolder+'_rep'+str(rep)
    print(f"Starting equilibration for: {output_subfolder}")
    # Execute eq1 through eq5 in series
    for i in range(1, 6):
        all_done = False
        started = False
        equilibration_counter = 0
        while not all_done:

            output_file = os.path.join(output_subfolder, f"eq{i}.log")
            cmd = f'python running_file.py {output_subfolder} eq{i} eq{i} {os.path.dirname(path_to_q)}'

            if not started:
                #check if the log file exists and is completed
                status = log_file_reader(output_file)
                if status == 1:
                    print(f"eq{i} already completed for: {output_subfolder}")
                    all_done = True
                else:
                #start the run
                    print(f"Running command: {cmd}")
                    os.popen(cmd + ' &').read() # Run in background
                    time.sleep(1)  # Brief pause to ensure command starts properly
                    started = True
            else:
                #check if the log file exists and is completed
                status = log_file_reader(output_file)
                if status == 1:
                    print(f"eq{i} completed for: {output_subfolder}")
                    all_done = True
                elif status == 2:
                    print(f"eq{i} encountered an error for: {output_subfolder}. Restarting...")
                    os.popen(cmd + ' &').read() # Run in background
                    equilibration_counter = 0
                else:
                    print(f"eq{i} still running for: {output_subfolder}. Waiting...")
                    equilibration_counter += 1
                    if equilibration_counter >= eq_wait_timer*60:  # After some time, die
                        print(f"eq{i} is taking too long for: {output_subfolder}. Exiting...")
                        sys.exit(1)
                    time.sleep(10)  # Wait before checking again
        print(f'equilibration step {i} done for: {output_subfolder}')
    print(f"Equilibration completed for: {output_subfolder}, {sim_type}")
        


#now we define the generic eq6_0_n and eq6_1_n runs
def generic_eq6_function(subfolder, index_of_run, sim_type, rep):
    

    output_subfolder = './neq_'+sim_type +'/'+ subfolder+'_rep'+str(rep)
    print(f"Starting generic eq6 run {index_of_run} for: {output_subfolder}")
    # Execute eq6_0_n and eq6_1_n in series
    for i in range(2):  # 0 and 1
        #modify the default eq6_0 and eq6_1 files to have the correct steps, as well as the correct input
        input_file = os.path.join(output_subfolder, f"eq6_{i}_default.inp")
        output_file = os.path.join(output_subfolder, f"eq6_{i}_rep"+str(index_of_run)+".inp")
        if index_of_run == 0:
            replacements = {
                "final       ": f" eq6_{i}_rep0.re",
                "energy      ": f" eq6_{i}_rep0.log",
                "restart     ": " eq5.re",
            
                "steps       ": f" {number_of_steps_eq5}"
            }
            
        else:
            replacements = {
                "final       ": f" eq6_{i}_rep{index_of_run}.re",
                "energy      ": f" eq6_{i}_junk.log",
                "restart     ": f" eq6_{i}_rep{index_of_run-1}.re",
                
                "steps       ": f" {number_of_steps_eq6}"
            }
        copy_and_modify_file(input_file, output_file, replacements)
        print(f"Processed: {output_file}")

        #no need to check if completed, as that is handled in the index of run 
        #we just need to check if started, start it if not, and wait until completion
        all_done = False
        started = False
        generic_eq6_counter = 0
        while not all_done:
            output_file = os.path.join(output_subfolder, f"eq6_{i}_{index_of_run}.log")
            if not started:
                cmd = f'python running_file.py {output_subfolder} eq6_{i}_rep{index_of_run} eq6_{i}_{index_of_run} {os.path.dirname(path_to_q)}'
                print(f"Running command: {cmd}")
                os.popen(cmd + ' &').read() # Run in background
                time.sleep(1)  # Brief pause to ensure command starts properly
                started = True
            else:
                status = log_file_reader(output_file)
                if status == 1:
                    print(f"Generic eq6_{i} run {index_of_run} completed for: {output_subfolder}")
                    all_done = True
                elif status == 2:
                    print(f"Generic eq6_{i} run {index_of_run} encountered an error for: {output_subfolder}. Exiting...")
                    sys.exit(1)
                else:
                    print(f"Generic eq6_{i} run {index_of_run} still running for: {output_subfolder}. Waiting...")
                    generic_eq6_counter += 1
                    if generic_eq6_counter >= eq6_generic_wait_timer*60:  # After some time, die
                        print(f"Generic eq6_{i} run {index_of_run} is taking too long for: {output_subfolder}. Exiting...")
                        sys.exit(1)
                    time.sleep(10)  # Wait before checking again
    print(f"Generic eq6 run {index_of_run} completed for: {output_subfolder}")

#and now the generic neq_0_n and neq_1_n runs
def neq_function(subfolder, index_of_run, sim_type, rep):
   

    output_subfolder = './neq_'+sim_type +'/'+ subfolder+'_rep'+str(rep)
    print(f"Starting neq run {index_of_run} for: {output_subfolder}")
    # Execute neq_0_n and neq_1_n in series
    for i in range(2):  # 0 and 1
        #modify the  files to have the correct steps, as well as the correct input
        input_file = os.path.join(output_subfolder, f"eq6_{i}_default.inp")
        output_file = os.path.join(output_subfolder, f"neq_{i}_rep"+str(index_of_run)+".inp")
        replacements = {
            "steps             ": f" {number_of_steps_neq}",
            "final       ": f" neq_{i}_junk.re",
            "energy      ": f" neq_{i}_junk.en",
            "restart     ": f" eq6_{i}_rep{index_of_run}.re",
            
            
        }
        copy_and_modify_file(input_file, output_file, replacements)
        append_lines_to_file(output_file, lines_to_add)

        print(f"Processed: {output_file}")
        #no need to check if completed, as that is handled in the index of run
        #we just need to check if started, start it if not, and wait until completion
        all_done = False
        started = False
        neq_counter = 0 
        while not all_done:
            output_file = os.path.join(output_subfolder, f"neq_{i}_{index_of_run}.log")
            if not started:
                cmd = f'python running_file.py {output_subfolder} neq_{i}_rep{index_of_run} neq_{i}_{index_of_run} {os.path.dirname(path_to_q)}'
                print(f"Running command: {cmd}")
                os.popen(cmd + ' &').read() # Run in background
                time.sleep(1)  # Brief pause to ensure command starts properly
                started = True
            else:
                status = log_file_reader(output_file)
                if status == 1:
                    print(f"neq_{i} run {index_of_run} completed for: {output_subfolder}")
                    all_done = True
                elif status == 2:
                    print(f"neq_{i} run {index_of_run} encountered an error for: {output_subfolder}")
                    all_done = True
                else:
                    print(f"neq_{i} run {index_of_run} still running for: {output_subfolder}. Waiting...")
                    neq_counter += 1
                    if neq_counter >= neq_wait_timer*60:  # After some time, die
                        print(f"neq_{i} run {index_of_run} is taking too long for: {output_subfolder}. Exiting...")
                        sys.exit(1)
                    time.sleep(10)  # Wait before checking again
    print(f"neq run {index_of_run} completed for: {output_subfolder}")

#this is the main function that handles the logic
#checks if the equilibration has been completed, then checks for eq6
#then it checks for the first incomplete eq_6_n run and starts the eq6_n - neq_n runs until the number of reps is reached
def runner_function(subfolder, sim_type, rep):
    
    
    script_prep_function(subfolder, sim_type, rep)
    equilibration_function(subfolder, sim_type, rep)
    #now we do the generic eq6 and neq runs. First, we identify the last file
    for index_of_run in range(0,number_of_reps):
        all_eq6_neq_done = True
        output_subfolder = './neq_'+sim_type + '/'+subfolder+'_rep'+str(rep)
        neq0_log = os.path.join(output_subfolder, f"neq_0_{index_of_run}.log")
        neq1_log = os.path.join(output_subfolder, f"neq_1_{index_of_run}.log")
        status_neq0 = os.path.exists(neq0_log)
        status_neq1 = os.path.exists(neq1_log)
        print(neq0_log)
        if not(status_neq0 or status_neq1):
            all_eq6_neq_done = False
            print(f'final index identified {index_of_run}')
            break
    if all_eq6_neq_done:
        print(f"All neq runs for index {index_of_run} are already completed. Moving to next index.")
        
    else:
        print(f"Starting generic eq6 and neq runs for index {index_of_run}.")
        for index_of_run in range(index_of_run, number_of_reps):
            if index_of_run % 10 == 0 and index_of_run > 20:
                clean_up_function(subfolder, sim_type, index_of_run, rep)
            
            generic_eq6_function(subfolder, index_of_run, sim_type, rep)
            neq_function(subfolder, index_of_run, sim_type, rep)




if __name__ == "__main__":
    if len(sys.argv) != 4:
        print("Usage: python folder_execution_neq.py <subfolder> <sim_type> <rep>")
        sys.exit(1)
    subfolder = sys.argv[1]
    sim_type = sys.argv[2]  # "water" or "protein"
    rep = int(sys.argv[3])  # the rep number
    print(f'Processing subfolder: {subfolder}, type {sim_type}, rep {rep}')
    runner_function(subfolder, sim_type, rep)
    print("All tasks completed.")


