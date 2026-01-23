#!/bin/bash
#
# Run CUDA implementation (qgpu) for vdW arithmetic rule validation
#
# This script uses the Python preprocessing pipeline (src/bin/qdyn.py) to:
# 1. Convert .top/.inp/.fep files to CSV format
# 2. Run the CUDA binary (qdyn_main) with the CSV files
#
# NOTE: The CUDA binary requires Linux with NVIDIA GPU. On macOS, use
# --preprocess-only to validate the CSV generation.
#
# Usage: ./run_cuda.sh [--preprocess-only]

set -e

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
Q_ROOT="$(cd "$SCRIPT_DIR/../../.." && pwd)"

QDYN_PY="$Q_ROOT/src/bin/qdyn.py"
QDYN_MAIN="$Q_ROOT/src/bin/qdyn_main"
WORKDIR="cuda_workdir"

PREPROCESS_ONLY=false
if [ "$1" = "--preprocess-only" ]; then
    PREPROCESS_ONLY=true
fi

if [ ! -f "$QDYN_PY" ]; then
    echo "Error: Cannot find $QDYN_PY"
    exit 1
fi

# Check binary architecture
check_binary_compat() {
    if [ ! -f "$QDYN_MAIN" ]; then
        echo "Warning: CUDA binary not found at $QDYN_MAIN"
        return 1
    fi

    local file_type=$(file "$QDYN_MAIN")
    local os_type=$(uname -s)

    if [[ "$os_type" == "Darwin" ]] && [[ "$file_type" == *"ELF"* ]]; then
        echo "Warning: CUDA binary is Linux ELF format, but running on macOS"
        echo "The binary cannot execute on this system."
        return 1
    fi

    if [ ! -x "$QDYN_MAIN" ]; then
        echo "Warning: CUDA binary is not executable"
        return 1
    fi

    return 0
}

cd "$SCRIPT_DIR"

echo "QGPU Arithmetic vdW Validation Test"
echo "===================================="
echo "Topology: dualtop.top"
echo "MD input: eq1.inp"
echo "FEP file: FEP1.fep"
echo "Workdir:  ./$WORKDIR"
echo ""

# Check if we can run the binary
CAN_RUN_BINARY=true
if ! check_binary_compat; then
    CAN_RUN_BINARY=false
    echo ""
fi

if [ "$PREPROCESS_ONLY" = true ]; then
    echo "Running preprocessing only (--preprocess-only flag set)..."
    echo ""
elif [ "$CAN_RUN_BINARY" = false ]; then
    echo "Continuing with preprocessing only..."
    echo "To run the full test, use a Linux system with NVIDIA GPU."
    echo ""
fi

# Run the preprocessing pipeline
# Note: The qdyn.py script will fail when trying to execute the binary on macOS,
# but the preprocessing (CSV generation) completes before that step.
set +e  # Don't exit on error - we want to check if preprocessing worked
python3 "$QDYN_PY" \
    --top dualtop.top \
    --md eq1.inp \
    --fep FEP1.fep \
    --workdir "$WORKDIR" \
    --core q7-gpu \
    --verbose 2>&1 | tee eq1_cuda.log
QDYN_EXIT=$?
set -e

# Check if preprocessing succeeded by looking for topo.csv
if [ -f "$WORKDIR/dualtop/topo.csv" ]; then
    echo ""
    echo "Preprocessing completed successfully."
    echo "CSV files generated in: $WORKDIR/dualtop/"

    # Verify vdw_rule is set to 2 (arithmetic)
    VDW_RULE=$(sed -n '9p' "$WORKDIR/dualtop/topo.csv")
    if [ "$VDW_RULE" = "2" ]; then
        echo "PASS: vdw_rule = 2 (arithmetic) correctly set in topo.csv"
    else
        echo "FAIL: vdw_rule = $VDW_RULE (expected 2)"
        exit 1
    fi
else
    echo ""
    echo "Error: Preprocessing failed - topo.csv not generated"
    echo "Check eq1_cuda.log for details."
    exit 1
fi

# Check if full run completed
if [ "$PREPROCESS_ONLY" = true ] || [ "$CAN_RUN_BINARY" = false ]; then
    echo ""
    echo "Preprocessing validation complete."
    echo "To run the full CUDA simulation, execute on Linux with NVIDIA GPU."
    exit 0
fi

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
