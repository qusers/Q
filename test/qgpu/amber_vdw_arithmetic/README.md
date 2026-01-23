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
| `run_fortran.sh` | Run qdyn6 (Fortran) |
| `run_cuda.sh` | Run qgpu (CUDA) |
| `compare_energies.sh` | Compare outputs |
| `compare_energies.py` | Python energy comparison script |

## Running the Test

```bash
cd test/qgpu/amber_vdw_arithmetic

# 1. Run Fortran reference
./run_fortran.sh ../../../bin/qdyn6

# 2. Run CUDA implementation
./run_cuda.sh ../../../bin/qgpu

# 3. Compare and validate
./compare_energies.sh
```

## Validation Criteria

Step 0 energies must match within **0.01 kcal/mol**:

- Total potential energy
- vdW energy
- Coulomb energy
- Q-surr energies (vdW and electrostatic)

## Expected Output

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
