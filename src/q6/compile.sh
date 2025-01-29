#!/bin/bash

if [ $# -ne 1 ]; then # Check if an argument was provided
    echo "Error: Please provide a compiler configuration"
    echo "Usage: ./compile_q.sh <compiler>"
    echo "Available compiler options: gcc, epyc4, ifort, osx, pgi"
    exit 1
fi

COMPILER=$1

echo "Compiling qprep & qfep with default COMP=gcc configuration..."
make qprep qfep COMP=$COMPILER

if [ $? -ne 0 ]; then # check if command was successful
    echo "Error: Failed to compile qprep and qfep"
    exit 1
fi

echo "Compiling simulation modules (qdyn, qdum, qcalc) with COMP=$COMPILER..."
make qdyn qdum qcalc COMP=$COMPILER

if [ $? -ne 0 ]; then # check if command was successful
    echo "Error: Failed to compile qdyn, qdum, and qcalc"
    exit 1
fi

echo "Compiling MPI version with COMP=$COMPILER..."
make mpi COMP=$COMPILER

if [ $? -ne 0 ]; then # check if command was successful
    echo "Error: Failed to compile MPI version"
    exit 1
fi

echo "Compilation completed successfully!"