#!/usr/bin/env python3
"""
Compare energy outputs between Fortran (qdyn6) and CUDA (qgpu) implementations.

Extracts step 0 energies from both log files and validates they match within
the specified tolerance. This validates the arithmetic vdW combination rule
implementation in the CUDA kernel against the Fortran reference.
"""

import argparse
import re
import sys
from dataclasses import dataclass
from pathlib import Path


@dataclass
class EnergyValues:
    """Container for energy values extracted from a log file."""
    total: float = 0.0
    potential: float = 0.0
    kinetic: float = 0.0
    q_surr_el: float = 0.0
    q_surr_vdw: float = 0.0


def parse_log_file(log_path: Path) -> EnergyValues:
    """
    Parse a qdyn/qgpu log file and extract step 0 energies.

    Returns:
        EnergyValues with total, potential, kinetic, Q-surr el, and Q-surr vdw.
    """
    content = log_path.read_text()
    energies = EnergyValues()

    # Find the "Energy summary at step 0" block
    # Pattern matches lines between "Energy summary at step      0" and "==="
    energy_block_match = re.search(
        r'Energy summary at step\s+0.*?={3,}',
        content,
        re.DOTALL
    )

    if energy_block_match:
        block = energy_block_match.group()

        # Extract SUM line: "SUM      <total>    <potential>    <kinetic>"
        sum_match = re.search(r'SUM\s+([\d.-]+)\s+([\d.-]+)\s+([\d.-]+)', block)
        if sum_match:
            energies.total = float(sum_match.group(1))
            energies.potential = float(sum_match.group(2))
            energies.kinetic = float(sum_match.group(3))

    # Find Q-atom energies at step 0
    q_block_match = re.search(
        r'Q-atom energies at step\s+0.*?={3,}',
        content,
        re.DOTALL
    )

    if q_block_match:
        q_block = q_block_match.group()

        # Extract Q-surr line: "Q-surr. <state> <lambda> <el> <vdw>"
        # Take the first state (usually state 1)
        qsurr_match = re.search(
            r'Q-surr\.\s+\d+\s+[\d.]+\s+([\d.-]+)\s+([\d.-]+)',
            q_block
        )
        if qsurr_match:
            energies.q_surr_el = float(qsurr_match.group(1))
            energies.q_surr_vdw = float(qsurr_match.group(2))

    return energies


def compare_energies(
    fortran: EnergyValues,
    cuda: EnergyValues,
    tolerance: float = 0.01
) -> tuple[bool, list[str]]:
    """
    Compare energy values between Fortran and CUDA implementations.

    Returns:
        Tuple of (all_passed, list of difference messages)
    """
    comparisons = [
        ("Total energy", fortran.total, cuda.total),
        ("Potential energy", fortran.potential, cuda.potential),
        ("Kinetic energy", fortran.kinetic, cuda.kinetic),
        ("Q-surr electrostatic", fortran.q_surr_el, cuda.q_surr_el),
        ("Q-surr vdW", fortran.q_surr_vdw, cuda.q_surr_vdw),
    ]

    all_passed = True
    messages = []

    for name, f_val, c_val in comparisons:
        diff = abs(f_val - c_val)
        passed = diff <= tolerance
        status = "PASS" if passed else "FAIL"

        msg = f"{status}: {name:25s} Fortran={f_val:12.4f}  CUDA={c_val:12.4f}  Diff={diff:.6f}"
        messages.append(msg)

        if not passed:
            all_passed = False

    return all_passed, messages


def main():
    parser = argparse.ArgumentParser(
        description="Compare Fortran and CUDA energy outputs for vdW arithmetic rule validation"
    )
    parser.add_argument(
        "fortran_log",
        type=Path,
        help="Path to Fortran (qdyn6) log file"
    )
    parser.add_argument(
        "cuda_log",
        type=Path,
        help="Path to CUDA (qgpu) log file"
    )
    parser.add_argument(
        "--tolerance",
        type=float,
        default=0.01,
        help="Maximum allowed difference in kcal/mol (default: 0.01)"
    )

    args = parser.parse_args()

    # Check files exist
    if not args.fortran_log.exists():
        print(f"Error: Fortran log file not found: {args.fortran_log}")
        sys.exit(1)

    if not args.cuda_log.exists():
        print(f"Error: CUDA log file not found: {args.cuda_log}")
        sys.exit(1)

    # Parse log files
    print(f"Parsing Fortran log: {args.fortran_log}")
    fortran_energies = parse_log_file(args.fortran_log)

    print(f"Parsing CUDA log: {args.cuda_log}")
    cuda_energies = parse_log_file(args.cuda_log)

    # Compare energies
    print(f"\nComparing step 0 energies (tolerance: {args.tolerance} kcal/mol)")
    print("=" * 80)

    all_passed, messages = compare_energies(
        fortran_energies,
        cuda_energies,
        args.tolerance
    )

    for msg in messages:
        print(msg)

    print("=" * 80)

    if all_passed:
        print("\nRESULT: PASS - All energies match within tolerance")
        sys.exit(0)
    else:
        print("\nRESULT: FAIL - Some energies differ beyond tolerance")
        sys.exit(1)


if __name__ == "__main__":
    main()
