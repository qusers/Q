#!/bin/bash
#
# Compare Fortran and CUDA energy outputs
#
# Usage: ./compare_energies.sh
#
# Requires: eq1_fortran.log and eq1_cuda.log to exist (run run_fortran.sh and run_cuda.sh first)

set -e

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
cd "$SCRIPT_DIR"

FORTRAN_LOG="eq1_fortran.log"
CUDA_LOG="eq1_cuda.log"

# Check that both log files exist
if [ ! -f "$FORTRAN_LOG" ]; then
    echo "Error: $FORTRAN_LOG not found"
    echo "Run ./run_fortran.sh first"
    exit 1
fi

if [ ! -f "$CUDA_LOG" ]; then
    echo "Error: $CUDA_LOG not found"
    echo "Run ./run_cuda.sh first"
    exit 1
fi

# Run the Python comparison script
python3 compare_energies.py "$FORTRAN_LOG" "$CUDA_LOG" "$@"
