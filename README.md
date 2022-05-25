# Q-GPU #
Version control of **Q-GPU**, an adaptation of **Q** version 5.06 running on GPUs.

**Q** is a set of Molecular Dynamics (MD) tools tailored to the following specific kinds of free energy calculations:

1. Free Energy Perturbation (FEP)
2. Empirical Valence Bond (EVB)
3. Linear Interaction Energies (LIE)

This version includes a translation of the original **Q** fortran code to C/CUDA and Python.


## Authors ##
Chiel Jespers, Willem Jespers, Mauricio Esguerra, Johan Åqvist, Hugo Gutiérrez‐de‐Terán


## Installation ##
The frontend is built on Python and will run in versions > 3.6. It mainly uses native python libraries and only needs numpy as additional package with no further dependencies.

To compile the qdyn engine source code, you need a CUDA compiler. The code has been tested with the following versions:

- CUDA/10.1.243

To succesfully install and compile the code:

```bash
unset SSH_ASKPASS
mkdir ~/software
cd ~/software
git clone https://yourgitusernamehere@github.com/qusers/qgpu.git 
cd qgpu/src/core
make
```

After succesful compilation of **Q-GPU** you have to add the program to your system path by modifying your shell initiation script. 
If your shell is bash, you can add the following lines to your .bashrc file using a text editor. The following assumes that your user name is "johndoe" and the home directory is "/Users/johndoe/":

```bash
SOFT=/Users/johndoe/software
export QDIR=$SOFT/qgpu
export PATH=$QDIR/bin:$QDIR/src:$PATH  
```
Where $SOFT will be the place where your software folder is located at, e.g. /Users/johndoe/software

Once the q binaries are declared in your path you should be able to call all q binaries from your terminal.
To test that the path to your compiled **Q** binaries has been correctly assigned you can issue the following commands in the terminal:

```bash
source ~/.bashrc
env | grep qgpu

QDIR=/Users/johndoe/software/qgpu
```

Additiontally you can search for the main **Q-GPU** binary file with:

```bash
which qdyn
```


## NOTE to the current version ##
The Qprep tool from **Q** is needed for the preparation of molecular topology files required by the MD engine Qdyn. Currently, Qprep is provided as fortran code, which compiles on CPUs. The workflow for a **Q-GPU** free energy simulation consists then in:

- An initial topology preparation stage that runs on a regular CPU  
- MD sampling using Qdyn, which runs on a CUDA-based GPU  
- The FEP analysis tool (qfep) provided in python (running both in GPU or CPU)  


## Troubleshooting ##
If you receive error messages during compilation please report them to the program authors including the compiler used (e.g. CUDA), the compiler version (e.e. 10.1.243), and the error message.


## Testing ##
**Q-GPU** includes various tests that compare the output of the original fortran code with the C/CUDA code. They are situated in the test folder and include:

1. interactions  
2. physical-properties  

The first folder includes test cases for the different type of interactions in **Q**, that is water-water (w-w), solute-solute (p-p) and Qatom-Qatom (q-q) interactions, and any mixture thereof.
These tests run a single point energy calculation and are compared with the output from Q5.07. The tests can be run separately following the instructions in each folder, or all at once using the run_test.py script (TODO!).

In the second folder, we provide test cases for the calculation of solvation free energies of side-chain mimics, and several protein-ligand binding cases (CDk2 and A2aAR, TODO!). The details for such calculations are described in our QligFEP paper:

- Jespers et al. (<https://doi.org/10.1186/s13321-019-0348-5>).

## Benchmarking ##

We have included a benchmark set of water spheres of sizes 10-30A (in increments of 5). Table generated with https://www.tablesgenerator.com/markdown_tables

| sphere | cpu Intel(R) Xeon(R) CPU E5-2650 v4 @ 2.20GHz Time in seconds | gpu NVIDIA GeForce GTX 1080  Time in seconds |
|--------|:-------------------------------------------------------------:|:--------------------------------------------:|
| 10A    |                                                         6.838 |                                        1.988 |
| 15A    |                                                        60.698 |                                        4.882 |
| 20A    |                                                       368.657 |                                       19.233 |
| 25A    |                                                      1257.150 |                                       59.948 |
| 30A    |                                                      4060.083 |                                      192.180 |

# VERSION NOTES: #

**19/08/2020**  
Generating first version of **Q-GPU** readme.  
