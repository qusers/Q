#!/bin/bash
#
# Run Fortran reference (qdyn6) for vdW arithmetic rule validation
#
# Usage: ./run_fortran.sh /path/to/qdyn6

set -e

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"

if [ -z "$1" ]; then
    echo "Usage: $0 /path/to/qdyn6"
    echo ""
    echo "This script runs the Fortran (qdyn6) implementation to generate"
    echo "reference energies for validating the CUDA implementation."
    exit 1
fi

QDYN6="$1"

if [ ! -x "$QDYN6" ]; then
    echo "Error: Cannot execute $QDYN6"
    echo "Make sure the path is correct and the binary has execute permission."
    exit 1
fi

cd "$SCRIPT_DIR"

echo "Running Fortran (qdyn6) on eq1.inp..."
"$QDYN6" eq1.inp > eq1_fortran.log

if grep -q 'terminated normally' eq1_fortran.log; then
    echo "Fortran run completed successfully."
    echo "Output saved to: eq1_fortran.log"
else
    echo "Error: Fortran run did not terminate normally."
    echo "Check eq1_fortran.log for details."
    exit 1
fi
