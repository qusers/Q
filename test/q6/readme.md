Description of Q tests
======================

## Test 1 -- Ligand in water plain MD. ##

This ligand has pubchem cid 2203539:  
<https://pubchem.ncbi.nlm.nih.gov/compound/2203539>

And its canonical SMILES is: 
C1=CC(=CN=C1)CN2C(=O)C3=C(C2=O)C=C(C=C3)C(=O)NC4=CC(=CC(=C4)C(=O)O)C(=O)O


## Test 2 -- Same as test 1 but smoother ##

This test is the  same as in test 2, that is,  the 2203539 ligand, but
in this case a smoother  equilibration protocol is used where heating,
relaxation, and final equilibration steps are decoupled.


## Test 3 -- A general FEP example ##

The idea is to have a FEP test here. 
Nothing here yet, but, perhaps Johan's sodium could go here.


## Test 4 -- A check for a soft-pair failure ##

An error was reported by Matt Mills on using soft-pairs.
Here we try to reproduce the error.


## Test 5 -- DiHydroFolateReductase (DHFR) ##

The CHARMM-AMBER benchmark molecule.  DHFR has now become the standard
benchmark for  explicit solvent molecular dynamics in  the majority of
equations  of motion  solvers.  As  of  March 2018  John Chodera  from
OpenMM claims that  506 ns/day are done in this  test, and Ross Walker
from  AMBER claims 517.12  ns/day. The  big difference  is in  the GPU
cards they use. Chodera says the benchmark is using a GTX-1080Ti ($699
list price) and Ross Walker reports using a Titan V100 $8000.


## Test 6 -- DNA and RNA using CHARMM36 ##

Simple MD  in Spherical Boundary  Conditions to test the  stability of
DNA and RNA with respect to Arnott's fiber models.


