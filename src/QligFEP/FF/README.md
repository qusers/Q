# On adding new parameters to the force field;

The forcefields provided in this repo are provided as `.lib` and `.prm` files. These files are used by Q to define the force field parameters. The `.lib` file contains the atom types and the `.prm` file contains the parameters for the force field.

## Charmm & OPLS ForceFields

To add extra parameters to the forcefield, you will need to make some conversion. For the `OPLS2015` and the `CHARMM36` forcefields, you can obtain the values required by the `.prm` file by using the function:

```python
import math

def vdw_calc(sig, eps):
    # calculate LJ types
    sig = float(sig)
    eps = float(eps)
    Avdw1 = math.sqrt(4*eps*(sig**12))
    Avdw2 = math.sqrt(4*eps*(sig**12))
    Bvdw1 = math.sqrt(4*eps*(sig**6))
    Avdw3 = (math.sqrt(4*eps*(sig**12)))/math.sqrt(2)
    Bvdw23 = (math.sqrt(4*eps*(sig**6)))/math.sqrt(2)

    # returns unparsed floats, consider changing
    return[Avdw1, Avdw2, Bvdw1, Avdw3, Bvdw23]
```

## Amber ForceField

However, if you're using the `AMBER14sb` forcefield, the parameters to be used in the `.prm` should be different, as they follow the header:

`R*1`, `R*2`, `epsilon1`, `R*3`, `epsilon2&3`, `mass`

In this case, follow this definition:
- `R*1` is the `sigma` value
- `R*2` is 0
- `epsilon1` is the `eps` value
- `R*3` is the same `sigma` value
- `epsilon2&3` corresponds to `sigma/2`
- `mass` is the mass of the atom.

The following Python function outputs the parameters in the format required by the `AMBER14sb.prm` file (just make sure you have matching indentation with the other parameters in the file):

```python
def array_calc(atom, sig, eps, mass = None):
    # calculate LJ types
    sig = float(sig)
    eps = float(eps)
    if mass is None:
        mass = 0.0
    return (
        f"{atom:<13}{sig:>10}{0.0:>11}"
        f"{eps:>11}{sig:>11}{eps/2:>11.6f}{mass:>11.2f}"
    )
```

## Adding extra ion parameters

If you're looking into extra parameters for ions, for example, you can use the parameters from `TIP3P water`. For this, you should first [install Amber](https://ambermd.org/GetAmber.php#ambertools), so that you can use the `tleap` program. Then, you can load the `TIP3P` water model and investigate the parameters for your ions of interest. For example, to load the `TIP3P` water model and investigate the parameters for the `Na+` ion, you can use the following commands:

```bash
tleap
source leaprc.water.tip3p
```

Where you will see the output:
``` text
----- Source: env/AmberTools23/dat/leap/cmd/leaprc.water.tip3p
----- Source of env/AmberTools23/dat/leap/cmd/leaprc.water.tip3p done
Loading library: env/AmberTools23/dat/leap/lib/atomic_ions.lib
Loading library: env/AmberTools23/dat/leap/lib/solvents.lib
Loading parameters: env/AmberTools23/dat/leap/parm/frcmod.tip3p
Reading force field modification type file (frcmod)
Reading title:
This is the additional/replacement parameter set for TIP3P water
Loading parameters: env/AmberTools23/dat/leap/parm/frcmod.ions1lm_126_tip3p
Reading force field modification type file (frcmod)
Reading title:
Li/Merz ion parameters of monovalent ions for TIP3P water model (12-6 normal usage set)
Loading parameters: env/AmberTools23/dat/leap/parm/frcmod.ionsjc_tip3p
Reading force field modification type file (frcmod)
Reading title:
Monovalent ion parameters for Ewald and TIP3P water from Joung & Cheatham JPCB (2008)
Loading parameters: env/AmberTools23/dat/leap/parm/frcmod.ions234lm_126_tip3p
Reading force field modification type file (frcmod)
Reading title:
Li/Merz ion parameters of divalent to tetravalent ions for TIP3P water model (12-6 normal usage set)
```

Then, you can check the content of your file of interest, e.g.: `ions1lm_126_tip3p` (see the true path in your system). There, you will see something along the lines;
``` text
Li/Merz ion parameters of monovalent ions for TIP3P water model (12-6 normal usage set)
MASS
Na+  22.99

NONBON
Na+   1.475  0.03171494    HFE set for TIP3P water from Li et al., JCTC, 2015, 11, 1645
```
The parameters for the `Na+` ion are `1.475` and `0.03171494` (in the order `sigma` and `eps`). You can use these parameters to add the `Na+` ion to the forcefield (in the case of the `AMBER14sb` forcefield, implemented as `SOD`).

For the reasons why such parameters are differently configured in the `prm` files Q, please refer to the [user's manual](http://qdyn.no-ip.org/documents/qman.pdf). Specifically, the following section (page 31):

The Lennard-Jones potential can be written either as $ \frac{A_{ij}}{r^{12}} - \frac{B_{ij}}{r^6}$ or $ \varepsilon_{ij} \left( \left( \frac{R^*_{ij}}{r_{ij}} \right)^{12} - 2 \cdot \left( \frac{R^*_{ij}}{r_{ij}} \right)^6 \right) $, using the geometric or arithmetic rules, respectively, to combine parameters for pairs of atom types. Treatment of 1-4 interactions (LJ and electrostatic) is specific for each force field.

## Adding extra AMBER parameters
After installing our environment, you should be able to find AMBER parameters under the path `~/micromamba//envs/AmberTools23/dat/leap/parm`.

Upon inspecting these files, you'll see `frcmod` files. The parameters from those files are converted to Q lib/prm formats. Just a few operations need to be done.

In the example, we convert `frcmod.phosaa14SB`:
```text
Updated dihedral parameters for phosphorylated amino acids for use with ff14SB

MASS
OP      16.00   0.434   Monoanionic Phosphate oxygen (analog to O2)
...

BOND
HO-OQ   553.0   0.960   from HO-OH

ANGL
OP-P -OP        140.0   119.90  from O2-P -O2

DIHE
N -CX-CT-CG       1    0.152093       0   -1

IMPR
CA-CG-CA-HA     1.100  180   2

NONB
OP      1.7493  0.2100          modified acc. to FEP solvation energies for phosphate monoanions
```

Which are converted to `AMBER14sb.prm`:
```text
[atom_types]
OP               1.7493        0.0     0.2100     1.7493    0.87465       16.0
! first value is sigma(after NONB), second is 0.0, third is epsilon, fourth is sigma again, fifth is sigma/2, sixth is mass

[bonds]
HO           OQ               1106.0       0.96
! here the original has "HO-OQ   553.0   0.960   from HO-OH"
! we then multiply the first value by 2 and keep the second

[angles]
OP           P            OP                280.0      119.9
! Here the original has "OP-P -OP        140.0   119.90  from O2-P -O2"
! we multiply the first value by 2 and keep the second

[torsions]
N            CX           CT           CG             0.152093  -1.0        0.0   1.0
N            CX           CT           CG             0.172788  -2.0      180.0   1.0
N            CX           CT           CG             0.490342  -3.0        0.0   1.0
N            CX           CT           CG             0.045151   4.0      180.0   1.0
! This value is stated as `DIHE` in the original file. see:
! "N -CX-CT-CG       1    0.152093       0   -1"
! "N -CX-CT-CG       1    0.172788     180   -2"
! "N -CX-CT-CG       1    0.490342       0   -3"
! "N -CX-CT-CG       1    0.045151     180    4"
! The first value remains the same, then we have the periodicity (-1 to 4, in this case), phase shift (angles), and number of paths (always 1 in this case)

[impropers]
CA           CG           CA           HA                  1.1      180.0
! In the original we have: "CA-CG-CA-HA     1.100  180   2"
! The first value remains the same, then we have the angle. values are kept the same
```

# On the ForceField parameters:

The force field .prm files follow, (here in the example: AmberFF) the pattern:

```text
[options]
name                           Q-Amber14SB
type                           AMBER
vdw_rule                       arithmetic !vdW combination rule (geometric or arithmetic)
scale_14                       0.8333 ! electrostatic 1-4 scaling factor
switch_atoms                   off
improper_potential             periodic
improper_definition            explicit

[atom_types]
...

! Ligand vdW parameters
! End ligand vdW parameters

[bonds]
...

! Ligand bond parameters
! End ligand bond parameters

[angles]
...

! Ligand angle parameters
! End ligand angle parameters

[torsions]
...

! Ligand torsion parameters
! End ligand torsion parameters

[impropers]
...

! Ligand improper parameters
! End ligand improper parameters
```

This tutorial will show how to convert prm files from other force fields in the openMM format to the Q format. Let's take the (Amber14sb force field)[https://github.com/openmm/openmm/blob/8.1.1/wrappers/python/openmm/app/data/amber14/protein.ff14SB.xml] as an example.

⚠️ Note: the following instructions are only to be followed fo creating a forcefield with the `arithmetic` vdw rule (as in the Amber FF).

Source for the explanation of the OpenMM force field files: http://docs.openmm.org/7.6.0/userguide/application/05_creating_ffs.html

## atom_types:

For this you can follow the explanation on either [Charmm & OPLS ForceFields](#charmm--opls-forcefields) or [Amber ForceField](#amber-forcefield) sections.

To extract the correct values from the OpenMM force field file, you will need to make some conversions. There, `sigma` is reported in $nm$ and `eps` in $kJ/mol$. You can use the following Python function to convert these values to the units used by Q:

<!-- Not sure why multiplying 5.612 does the trick??? -->
```python
def atom_type_forces_from_openmm(sig, eps):
     # convert eps from kJ/mol to kcal/mol
     sig = round(sig * 5.612, 4)
     # convert sig to A
     eps = round(eps * 23.9005736 / 100, 3)
     return sig, eps
```
## bonds

For this, you will have to convert the `HarmonicBondForce`, but OpenMM uses spring forces reported in $kJ/mol/nm^2$, but Q uses these forces in $kcal/mol/A^2$ rounded up to 1 decimal case.

To convert these values you can use the Python function:

```python
def convert_kj_mol_nm2_to_kcal_mol_a2(value_kj_mol_nm2):
     conversion_factor = 23.9005736
     value_kcal_mol_a2 = round(value_kj_mol_nm2 * conversion_factor / 10000, 1)
     return value_kcal_mol_a2
```

For the length, you should use * 10 same values reported on the `HarmonicBondForce` present in the OpenMM force field file, but make sure to mu

## angles

According to the Q manual, the `[angles]` part of the file define the harmonic angle parameters. However, other transformations are needed to convert the values from the OpenMM force field file to the Q format.

Since OpenMM uses angles as $radians$ and spring constants in $kJ/mol/rad^2$, you will need to convert these values to the units used by Q. You can use the following Python function to convert these values:

```python
import math

def radians_to_degrees(radians):
    return radians * 180 / math.pi

# values for C, C, O angles
angle = 2.0943951023931953
Q_converted_angle = round(radians_to_degrees(angle), 1)

k = 669.44
Q_converted_k = convert_kj_mol_nm2_to_kcal_mol_a2(k) * 100
```
## torsions

According to Q's manual, this part is defined with the following columns:
1) atom type 1 or 0 or ? to match any atom type
2) atom type 2
3) atom type 3
4) atom type 4 or 0 or ? to match any atom type
5) force constant = 1/2 * barrier height (kcal * mol^-1)
6) periodicity (number of maxima per turn). Add a minus sign before to indicate that more components follow on subsequent lines, i. e. for a torsion potential with multiple components all but the last component should be entered with negative periodicity.
7) phase shift ($\delta/2$ define the location of first maximum) ($\degree$)
8) number of paths

## impropers

According to Q's manual, this part is defined with the following.

columns:
1) atom type 1 or 0 or ? to match any atom type
2) atom type 2
3) atom type 3
4) atom type 4 or 0 or ? to match any atom type
5) force constant = kcal / mol / rad^2
6) equilibrium angle

However, despite the force here being reported as $kcal/mol/rad^2$, it seems like multiplying this value by 10 does reproduce the correct values in the Q force field file.

See script `conversion.py` in this directory for a guide on how to convert the OpenMM force field files to the Q format.

# add section on how to add from extra forcefield:
