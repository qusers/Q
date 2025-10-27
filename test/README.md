===
We have implemented several tests:
name : p-p
topology: benzene-vacuum.top
sphere size: 20
description: simulation of benzene molecule in vacuum

name: q-p-benzene
topology: Na-benzene-vacuum.top
FEP file: FEP_benzene.fep
sphere size: 20
description: simulation of benzene molecule in vacuum, but it is considered a q atom, also Na atom is set to 0

name: q-p_Na
topology: Na-benzene-vacuum.top
FEP file: FEP_Na.fep
sphere size: 20
description: simulation of Na atom in vacuum, but it is considered a q atom, also benze molecule is set to 0

name: q-p-w_benzene
topology: Na-benzene-water.top
FEP file: FEP_benzene.fep
sphere size: 20
description: simulation of benzene molecule in water, but it is considered a q atom, also Na atom is set to 0

name: q-p-w_Na
topology: Na-benzene-water.top
FEP file: FEP_Na.fep
sphere size: 20
description: simulation of Na atom in water, but it is considered a q atom, also benzene molecule is set to 0

name: q-q
topology: benzene-vacuum.top
sphere size: 20
description: simulation of benzene atom in vacuum, but it is considered a q atom

name : w-p
topology: benzene-water.top
sphere size: 20
description: simulation of benzene molecule in water

name : w-q
topology: benzene-water.top
sphere size: 20
description: simulation of benzene molecule in water, but it is considered a q atom

name: w-w
topology: water.top
sphere size: 20
description: simulation of water sphere

name: boundary
topology: ala_wat.top
sphere size: 15
description: simulation of polypeptide, with part over the sphere boundary that gets fixed

name: polypeptide
topology: ala_wat.top
sphere size: 15
description: simulation of polypeptide in 15A water sphere

name: polypeptide25
topology: ala_wat.top
sphere size: 25
description: simulation of polypeptide in 25A water sphere

===
We will add the following tests:
name: thrombin
description: early equilibration of thrombin target (1a_1c)

name: thrombin_eq5
description: late equilibration of thrombin target

name: thrombin_0744_0256
description: mid perturbation of thrombin target

name: thrombin_0998_0002
description: late perturbation of thrombin target

Similar for other (harder?) targets

===
TODO for tests:
* A way to convert restart file to test scenario
* tee and tail command line options as described in the github issue