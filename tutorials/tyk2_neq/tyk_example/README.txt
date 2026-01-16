This directory contains scripts and configuration files for running
nonequilibrium molecular dynamics (NEQ MD) simulations to compute free
energy differences between two states (state 0 and state 1).

Typical condensed workflow: 
prepare the inputs using QligFEP v2, so they are ready for MD simulations
change the '1.protein' and '2.water' folder names to be 'protein' and 'water', and move them to be in the script folders
change the Q directory in folder_execution_neq.py if required, to point to you qdyn_neq executable
change the parameters in parameters.txt
sbatch j1.sh and j2.sh
wait until completion
run analyzer.sh
run combine_csv.py

Background
The workflow performs multiple independent nonequilibrium transformations
in both directions (0 → 1 and 1 → 0), records the work values, and analyzes
them using standard free energy estimators such as BAR.

The setup is intended for use with the Qdyn_neq / Q6  engine and is designed
to run on HPC systems using SLURM.

Instead of slowly transforming one state into another, NEQ MD performs
the transformation over a finite amount of time while recording the work
done on the system.

Each transformation produces a work value rather than a direct free energy.
By running many independent transformations in both directions (0 → 1 and
1 → 0), the forward and reverse work distributions can be combined to
estimate the free energy difference.

Running both directions is important because it improves statistical
efficiency and allows the use of BAR, which gives a robust estimate of the
free energy change.

The transformation is controlled by a time-dependent lambda parameter.
The rate at which lambda changes can be linear or sigmoidal, depending on
the chosen settings.


Parameter File

The parameter file defines how the NEQ transformations are carried out.
It is read by the main execution script and applies to all edges and
replicates.

Important parameters:

number_of_reps
Number of independent NEQ forward and reverse simulations per folder.

number_of_steps_eq5
Number of equilibration steps before the first NEQ simulation.
Increase this if the two states differ significantly.

number_of_steps_eq6
Number of equilibration steps between NEQ simulations.
Recommended to be greater than 250 steps.

number_of_steps_neq
Length of each NEQ simulation in MD steps.
Recommended to be at least 16000 steps.

L_parameter
Controls the steepness of the sigmoidal lambda schedule.
Larger values spend more time near lambda = 0 and lambda = 1.
Typical values are between 4 and 16.
If a linear lambda schedule is desired, set L_parameter to a very small
value (for example 0.0001) to avoid code changes.

Default values are used automatically if the parameter file is missing.

NEQ Execution Script, usage is typically though j1.sh (water) and j2.sh (protein)
CPU usage should be scaled roughly as:
(number of edges) × (number of replicates)
Threading is disabled for Python jobs by setting OMP_NUM_THREADS=1.

The underlying Python script runs the NEQ simulations for a given folder and
simulation type (protein or water).

The script:

Calls the Qdyn NEQ executable
Controls the time step
Optionally randomizes the random seed
Monitors job progress and terminates runs that exceed time limits

Important settings include:
Path to the Qdyn NEQ executable
MD time step
Wall-time limits for equilibration and NEQ stages


BAR analysis script - executed through analyzer.sh
The underlying python script does the following:

For a selected replicate number, reads work values from protein and water simulations
Performs BAR analysis
Writes per-replicate results to CSV files
Plots distributions of the work histograms and BAR uncertainty in the folder
python analysis_single_group.py <folder>
better to use analyzer.sh to call it in parallel though

Combining replicate outputs (combine_csv.py)

This script combines multiple replicate CSV files into a single CSV file.

Finds all files matching rep*.csv in the current directory
Uses the first column as the index
Extracts work values from the third column
Writes a combined output file named <prefix>_combined.csv
This is useful for aggregating replicate data for all the attempts


Work histogram and time evolution plotter:

This script plots cumulative work versus time for forward and reverse transformations
python analysis_histo_plot.py <log_folder>

These plots are useful for checking convergence and diagnosing poor
overlap between forward and reverse work distributions.

Archiving Script

A separate tarballing script is provided to compress NEQ output folders
after analysis in order to save disk space.

Typical Workflow

Edit the parameter file to set NEQ length, equilibration steps, andlambda schedule.
Submit SLURM jobs to run NEQ simulations for protein and water systems

example 'sbatch j1.sh' 'sbatch j2.sh'

Run BAR analysis and inspect work distributions.
Combine replicate CSV files as needed.
Archive NEQ folders once analysis is complete.

Notes and Practical Tips

Increase the NEQ length if forward and reverse work distributions show poor overlap.
Adjust the L parameter if the work gradients are steep.
Always inspect work distributions before trusting free energy estimates.
Balance computational cost by trading off number of replicates against NEQ length.