Q-GPU
=======
Version control of Q-GPU, an adaption of Q version 5.06 that runs on GPUs.

Q is a set of Molecular Dynamics (MD) tools tailored to the following specific kinds of free energy calculations:

1. Free Energy Perturbation (FEP)
2. Empirical Valence Bond (EVB)
3. Linear Interaction Energies (LIE)

This version includes a translation of the original Q fortran code to C/CUDA and Python.

## Authors:
Chiel Jespers, Willem Jespers, Johan Åqvist, Hugo Gutiérrez‐de‐Terán

## Installation
The Python environment works with Python > 3.6. It solely uses native python libraries, so no additional libraries will need to be installed.

To compile the qdyn engine source code, you need a CUDA compiler. The code has been tested with the following versions:

- CUDA/10.1.243

To succesfully install and compile the code:

```bash
unset SSH_ASKPASS
git clone https://yourgitusernamehere@github.com/qusers/qgpu.git
cd qgpu/src/core
make
```

After this you have to add the program to your system path by modifying your shell initiation script, that is, if your shell is bash, you can add the following lines to your .bashrc file using a text editor:

```bash
SOFT=/Users/johndoe/software
export QDIR=$SOFT/qgpu
export PATH=$QDIR/bin:$PATH  
```
Where $SOFT will be the place where your software folder is located at, e.g. /Users/johndoe/software

Once the q binaries are declared in your path you should be able to call all q binaries from your terminal.
To test that the path to your compiled Q binaries has been correctly assigned you can issue the following commands in the terminal:

```bash
source .bashrc
echo $path | grep qsource

/Users/johndoe/software/qsource
```

## Troubleshooting

If you receive error messages during compilation please report them to ???? including the compiler used (e.g. intel fortran), the compiler version, and the error message.

NOTES:
=========

19/08/2020

Generating first version of Q-GPU readme.
