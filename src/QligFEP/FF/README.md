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
def array_calc(sig, eps, mass):
    # calculate LJ types
    sig = float(sig)
    eps = float(eps)
    return f"{sig:.4f}, 0.0, {eps:.6f}, {sig:.4f}, {sig/2:.7f}, {mass:.3}"
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
