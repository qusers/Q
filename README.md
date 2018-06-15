Q for GPU development.
================================================================================

This a fork of Q6 with a large amount of differences from that
original code. This will form the starting point for GPU implementation using
OpenACC, though this might change to CUDA-fortran (Portland)
**Use at your own risk.**

The forked code is released under a GNU General Public License, Version 2.0
(GPLv2), that means that all modifications done here belong to the 
Public license owners, that is, Paul Bauer, j. marelius, j. aqvist,
s. c. l. kamerlin.

Status
------
[![Build Status](https://img.shields.io/travis/esguerra/Q6/master.svg?style=plastic)](https://travis-ci.org/esguerra/Q6)


The main idea of this fork is to aim for clean and readable Q code.  
The code is in the process of being reorganized and cleaned-up
according to best coding practices.  

We are starting to document the code using doxygen and ford.   
To generate the doxygen documentation do:  

    cd docs/developers  
    doxygen DoxygenConfigFortran  


To generate the ford docs:  
    
    cd docs/developers  
    ford -d ../../src -o devdocs forddocs.md  


Some basic development rules somewhat similar to those of the
**gromacs** developers team follow.


Coding Standards
----------------

The **GROMACS** developers have identified the following important main
points for taking into account when organizing a molecular dynamics
code.  

1.  **Code formatting** - how to indent code, how to start and end subroutines
    etc.  
2.  **Code constructs** - argument order, return values, encapsulation
    using abstract data types  
3.  **Interfaces** - the Application Programming Interface should say it
    all  
4.  **Comments in code** - comments in code that ford can use  
5.  **Compilation** - using different hardware  
6.  **Allowed Fortran Features**  
7.  **Error Handling**  
8.  **Benchmarking**  
9.  **Accuracy Testing**   

For fortran maybe slightly different rules will be needed. The fluidity project 
is a good example to draw from.

  
1. Code formatting
------------------

* No tabs, spaces only.  
* Two spaces for indentation of each level.  
* No more than 80 characters to allow for easy code visualization
  across editors and screens. Specially important for easy mobile
  device coding, reading.  
* No trailing whitespaces. They are not seen by default in most
  editors but still count as changes in git.  
* Column 36 for : on variable declaration.  


2. Code constructs
------------------
*to do*  


3. Interfaces
-------------
*to do*  


4. Comments in code
-------------------
*to do*  


5. Compilation
--------------

5.1  gfortran from gcc
----------------------
At the moment this fork includes a few features of the fortran
2008 standard. The fortran 2008 standard is the latest fortran
standard. Work is being done on the fortran 2015 standard but it is
not yet released. Many features of fortran 2008 are already included
in the GNU-Compiler-Collection gfortran.

- For compiling Q a new version of gcc is needed, for now we are
building with gcc-7.3 though we will probably change soon to gcc-8.
A script is added to this repository to install gcc with offloading
to NVIDIA GPU's.

- **Fortran 2008**(New Features)  
Coarrays for parallel computing.  
Possibility to use submodules, and submodules of submodules. Of help
for very large programs.  
*do concurrent*  

5.2 compilation in windows 10
-----------------------------

You will  have to install  mingw64 in order  to compile Q  in windows.
Once  mingw64 and  gcc-fortran  are  installed you  can  just use  the
standard compilation for gcc of Q.  We recommend using msys2, and then
the mingw64 shell.  You will need  git to clone the repository, so you
will need to  install git using *pacman* which is  the package manager
of mingw64. You will also need make and gcc-fortran.


    pacman -S pacman mingw-w64-x86_64-gcc-fortran
    pacman -S git make
   
Then just clone this fork:

    git clone https://github.com/esguerra/Q6.git

And compile with:

    cd src
    make all COMP=gcc

We have tested with gcc 6.3.0 and the compilation works well. For
compilation of the parallel *qdynp* you will need MPI in windows. We
have not tested Windows MPI yet, nor the openMPI port for cygwin.



6. Allowed Fortran Features
---------------------------
*to do*  


7. Error Handling
-----------------
*to do*  


8. Benchmarking
---------------
*to do*  


9. Accuracy Testing
-------------------
*to do*  



