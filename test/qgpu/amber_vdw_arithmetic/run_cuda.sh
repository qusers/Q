#!/bin/bash
#
# Run CUDA implementation (qgpu) for vdW arithmetic rule validation
#
# This script uses the Python preprocessing pipeline (src/bin/qdyn.py) to:
# 1. Convert .top/.inp/.fep files to CSV format
# 2. Run the CUDA binary (qdyn_main) with the CSV files
#
# Usage: ./run_cuda.sh

set -e

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
Q_ROOT="$(cd "$SCRIPT_DIR/../../.." && pwd)"

QDYN_PY="$Q_ROOT/src/bin/qdyn.py"
QDYN_MAIN="$Q_ROOT/src/bin/qdyn_main"
WORKDIR="cuda_workdir"

if [ ! -f "$QDYN_PY" ]; then
    echo "Error: Cannot find $QDYN_PY"
    exit 1
fi

if [ ! -x "$QDYN_MAIN" ]; then
    echo "Error: Cannot execute $QDYN_MAIN"
    echo "Make sure the CUDA binary has been built and has execute permission."
    exit 1
fi

cd "$SCRIPT_DIR"

echo "Running QGPU preprocessing and CUDA binary..."
echo "Topology: dualtop.top"
echo "MD input: eq1.inp"
echo "FEP file: FEP1.fep"
echo "Workdir:  ./$WORKDIR"
echo ""

python3 "$QDYN_PY" \
    --top dualtop.top \
    --md eq1.inp \
    --fep FEP1.fep \
    --workdir "$WORKDIR" \
    --core q7-gpu \
    --verbose 2>&1 | tee eq1_cuda.log

if grep -q 'terminated normally' eq1_cuda.log; then
    echo ""
    echo "CUDA run completed successfully."
    echo "Output saved to: eq1_cuda.log"
else
    echo ""
    echo "Error: CUDA run did not terminate normally."
    echo "Check eq1_cuda.log for details."
    exit 1
fi
