# vdW Arithmetic Rule Validation Test

This test validates the arithmetic vdW combination rule implementation in the CUDA kernel (qgpu) against the Fortran reference implementation (qdyn6).

## Purpose

AMBER force fields use the arithmetic combination rule for vdW interactions (Lorentz-Berthelot), while OPLS uses the geometric rule. This test ensures the CUDA implementation correctly handles the arithmetic rule (`vdw_rule = 2` / `FF_TYPE 2`) by comparing energies against the established Fortran code.

## System Description

- **Protein**: Tyk2 kinase
- **Ligands**: ejm_31 → ejm_42 FEP pair
- **Force field**: AMBER14sb
- **Atoms**: ~7600 total (protein + water + 2 ligands)
- **Sphere radius**: 25 Å

The topology file (`dualtop.top`) contains `FF_TYPE 2`, indicating arithmetic vdW combination.

## Files

| File | Description |
|------|-------------|
| `eq1.inp` | MD input file (10000 steps, fixed random seed) |
| `dualtop.top` | Topology with AMBER14sb + dual ligand |
| `FEP1.fep` | FEP definition (atom types, charges, softcore) |
| `AMBER14sb_ejm_31_ejm_42_merged.prm` | Merged parameter file |
| `ejm_31.lib` | Ligand 1 library |
| `ejm_42_renumber.lib` | Ligand 2 library (renumbered) |
| `run_fortran.sh` | Run qdyn (Fortran) |
| `run_cuda.sh` | Run qgpu (CUDA) via Python preprocessing |
| `compare_energies.sh` | Compare outputs |
| `compare_energies.py` | Python energy comparison script |

## Requirements

- **Fortran (qdyn)**: macOS or Linux, compiled from `src/q6/`
- **CUDA (qgpu)**: Linux with NVIDIA GPU and CUDA toolkit
  - The CUDA binary (`src/bin/qdyn_main`) is Linux ELF format
  - On macOS, only preprocessing validation is possible

## Running the Test

### On Linux (full test)

```bash
cd test/qgpu/amber_vdw_arithmetic

# 1. Run Fortran reference
./run_fortran.sh ../../../src/q6/bin/q6/qdyn

# 2. Run CUDA implementation (requires NVIDIA GPU)
./run_cuda.sh

# 3. Compare and validate
./compare_energies.sh
```

### On macOS (preprocessing validation only)

```bash
cd test/qgpu/amber_vdw_arithmetic

# 1. Run Fortran reference (works on macOS)
./run_fortran.sh ../../../src/q6/bin/q6/qdyn

# 2. Validate preprocessing generates correct CSV files
./run_cuda.sh --preprocess-only

# This verifies that vdw_rule=2 is correctly written to topo.csv
```

## Validation Criteria

Step 0 energies must match within **0.01 kcal/mol**:

- Total potential energy
- vdW energy
- Coulomb energy
- Q-surr energies (vdW and electrostatic)

## Expected Output

### Preprocessing validation (macOS)

```
QGPU Arithmetic vdW Validation Test
====================================
...
Preprocessing completed successfully.
CSV files generated in: cuda_workdir/dualtop/
PASS: vdw_rule = 2 (arithmetic) correctly set in topo.csv
```

### Energy comparison (Linux)

```
Comparing step 0 energies (tolerance: 0.01 kcal/mol)
================================================================================
PASS: Total energy               Fortran=   -1234.5678  CUDA=   -1234.5679  Diff=0.000100
PASS: Potential energy           Fortran=   -2345.6789  CUDA=   -2345.6790  Diff=0.000100
...
================================================================================

RESULT: PASS - All energies match within tolerance
```

## Troubleshooting

If the test fails:

1. Check that both binaries were compiled with the same precision (double recommended)
2. Verify `vdw_rule = 2` appears in both log files during parameter parsing
3. Look for NaN or infinity values in energy outputs
4. Check that the CUDA kernel is using the arithmetic combination formula:
   - σ_ij = (σ_i + σ_j) / 2
   - ε_ij = sqrt(ε_i * ε_j)

## Implementation Details

The vdW combination rule is controlled by:


**Fortran**: Right after section defining the number of atom types (`= No. of atom types`). For example:
- Geometric: `1 = vdW combination rule (1 = Geom. / 2 = Arit.)` (see [example](../../q6/test1/nc12/lig_w.top))
- Arithmetic: `2 = vdW combination rule (1 = Geom. / 2 = Arit.)` (see [example](./dualtop.top))

**CUDA**: Python preprocessing writes either `1` or `2` to line 9 of `topo.csv` (generated with [topology.py](../../../src/Qgpu/topology.py))

The arithmetic rule is parsed from topology sections labeled `R* normal:` and `epsilon normal:` (as opposed to `sqrt (Aii) normal:` and `sqrt (Bii) normal:` for geometric).
