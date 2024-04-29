# QligFEP  

This folder contains a step by step tutorial to setup and
analyze QligFEP calculations, as reported in Jespers et al.
(https://doi.org/10.1186/s13321-019-0348-5).  

The workflow consists of four main steps:  

1. Generate ligand parameters from `.sdf` files using OpenFF & use Lomap to generate the perturbation mapping;
2. Prepare the protein to run under spherical boundary conditions;
3. Generate all the FEP related input files.
4. Analyse the results.

⚠️ WARNING ⚠️ This tutorial is by no means complete, and it's currenly under development. If you have any suggestions/contributions/comments to make, that's all very welcome. There are a few commends on [issue 22](https://github.com/qusers/Q/issues/22), and if you have any comments please add it there so we can have a roadmap of things to fix/add.

# ligprep

Generating ligand parameters takes a while. Therefore we advice running it on the background. For example, to generate the parameters for the ligands in `Tyk2_ligands.sdf`, run:
```bash
nohup qparams -i Tyk2_ligands.sdf > Tyk2_qparams.log 2>&1 &
```
Running this command will output `.lib`, `.prm`, and `.pdb` files for each ligand in the `.sdf` file. Q uses these files as input for the relative free energy calculation.

For the sake of the tutorial, we have already generated the parameters for the ligands in `Tyk2_ligands.sdf` using the same command. You can find them in the `ligparams` folder.

Finally, generate the perturbation mapping using lomap:
```bash
qlomap -i Tyk2_ligands.sdf
```

Lomap natively requires a directory with separate `.sdf` files for each of the ligands to be used as input. Whenever `qlomap` is called to process a single `.sdf` file, our wrapper will create the directory *on the fly* for you and place the `lomap.json` inside the created directory. The directory will be named after the input `.sdf` file, without the `.sdf` extension.

# Water sphere

If you're working with a protein structure for which you don't have a water molecule, you can generate that using `qprep_prot`. You can check the arguments for that CLI through:

```bash
qprep_prot -h
```

You will see that for that you need to calculate the center of geometry for the water sphere. This can be done by calling `qcog` to a ligand file (either a `.sdf` or a `.pdb` file). If you're running this command for a `.sdf` file with multiple ligands, the script will calculate the center of geometry for all the ligands and then output the averave of all of them. It's up to you to choose which center of geometry to use.

As an example of `qprep_prot` usage, you can run, for the protein file `/Q/tutorials/Tyk2/setupFEP/amber`, the following command:

```
qprep_prot -cog -4.727 26.163 -30.542 -pdb protein.pdb -FF AMBER14sb -log DEBUG
```

# qprep debugging

**TODO** - input from other users is welcome

# setupFEP

## Intro to naming conventions

Due to the different atom naming conventions on the `.lib` and `.prm` implementations of each of the [forcefields used in Q](https://github.com/qusers/Q/tree/main/src/QligFEP/FF), you will have to make sure that the atoms in your protein file match the respective `.lib` and `.prm` files within your protein.pdb file.

For the sake of this tutorial, this atom renaming was already done. Please copy the `water.pdb` file found in this directory to either `amber` or `opls`, depending on the forcefield you want to use in your calculations.

TODO: make this more extensive...

## Starting with FEP calculations

Now that we have all the files generated on the [ligprep](#ligprep) step, we can proceed to the next step: setting up the FEP calculations.

At the setupFEP directory, you'll find `protein.pdb` and `water.pdb`, both of which are needed to run the FEP calculations. Additionally, you'll need to copy:

1. The generated `lomap.json` file from `ligprep/Tyk2_ligands/`
2. The `.lib`, `.prm`, and `.pdb` files from the `ligparams` folder to the `setupFEP` folder.

You can do so by running:

```bash
cp ligparams/* setupFEP/
cp ligprep/Tyk2_ligands/lomap.json setupFEP/
```

Now, because of the different forcefield [naming conventions](#intro-to-naming-conventions), we have to proceed with a different protein structure depending on the forcefield that you want to use. At the moment, both OPLS2015 and AMBER14sb are demonstrated. See below.

### Using the OPLS2015 forcefield

Please copy the file we have under `/Q/tutorials/Tyk2/setupFEP/opls/protein.pdb` to your `setupFEP` directory. E.g. by doing (if you're in the setupFEP directory):

```bash
cp opls/protein.pdb .
```
After doing this, you should be able to setup the calculations by running:

```bash
setupFEP -FF OPLS2015 -c TETRA -r 25 -l 0.5 -R 10 -w 100 -j lomap.json -S sigmoidal -clean .dcd
```

### Using the AMBER14sb forcefield

Please copy the file we have under `/Q/tutorials/Tyk2/setupFEP/amber/protein.pdb` to your `setupFEP` directory.  E.g. by doing (if you're in the setupFEP directory):

```bash
cp amber/protein.pdb .
```

After doing this, you should be able to setup the calculations by running:

```bash
setupFEP -FF AMBER14sb -c TETRA -r 25 -l 0.5 -R 10 -w 100 -j lomap.json -S sigmoidal -clean .dcd
```

### For converting protein files prepared in Maestro to AMBER14sb's naming convention

If you're going to use the AMBER force field, you'll likely run into issues where the atom names in your `protein.pdb` file don't match the atom names in the forcefield's `.lib` file.

If you need to correct this for your own protein structure, you can use the command below to rename the atoms in your `protein.pdb` file within `Tyk2/setupFEP/opls/protein.pdb` so you have the correct `AMBER` atom naming:
```bash
pdb2amber -i protein.pdb -o protein_amber.pdb
```
*Note*: setupFEP **always** uses a `protein.pdb` file as input for the FEP calculations. If you decide to use the AMBER forcefield, you should rename the `protein_amber.pdb` file to `protein.pdb`.

*Tip*: always keep a backup of your original `protein.pdb` file, in case you change your mind and decide on a different forcefield later on.

# Job submission

Now that you ran `setupFEP`, you should have both `1.water` and `2.water` directories

# END OF THE UPDATED PART OF THE TUTORIAL

**TODO**: Update tutorials so that we have a full version;

This folder contains pregenerated 3D-coordinates of the ligand.
Note that these 3D coordinates are not generated by setupFEP, and
you are responsible to get reliable coordinates. It is important
to keep in mind that FEP handles well small changes in the system,
but if the phase space is not sufficiently overlapping the results
become increasingly unreliable. The ligands in this folder are
all overlaying a similar core region, and small substituents on
are investigated.  

For AMBER parameters, please refer to qtools:  

https://zenodo.org/badge/latestdoi/80016679  

Since we have quite a few ligands (n = 16), we will generate the
files using a script from the script folder, e.g. run:  

    python $qligfep/scripts/generate_opls.py

This will write out the .lib, .prm and .pdb files needed for
stage 3. The .log files are generated by ffld_server.  

-- Note --  
To successfully run this script, a working version of Maestro is
needed, since the OPLS parameters are derived from ffld_server
output. See the setup section in $setupFEP/README.md  

-- Note --  
In order to successfully run generate_charmm.py:  
- Have cgenff installed and properly referenced in the script  
- Have a running version of openbabel  
- make sure that you have .mol2 files with partial charges added  

-- Note --  
Due to an implementation error in **Q**, 4 letter atom names are not
properly written out. In the case of particularly halogen atoms
you would have to double check and adjust the .pdb file, e.g.
in the case of 17.pdb in this tutorial, change:  

HETATM   25  Br25 LIG   900       3.933  25.020   9.894  

to:  

HETATM   25  Br25LIG   900       3.933  25.020   9.894  

# 2.protprep  
The second step is needed to prepare a protein for simulations
under spherical boundary conditions (SBC). Note that this script
does not attempt to assign protonation states, so one needs to run
a preparation step first (e.g.  Maestro's Protein Preparation
Wizard). Alternatively, protonation states can be assigned manually
by using the explicit names as provided in $setupFEP/FF/*.lib files.

The center of the sphere can be based on a residue in the protein,
or as in this case by explicitly providing the center of geometry
(COG) coordinates from a ligand. E.g. to get those used in this
example run:

    python $qligfep/scripts/COG.py ../1.ligprep/17.pdb

Which returns the coordinates 0.535:26.772:8.819, which can be
directly put in protprep.py:

    python $qligfep/protprep.py -p 1h1s_PrepWiz.pdb -r 22 -c 0.535:26.772:8.819 -w -P CSB

You can run $qligfep/protprep.py -h for a full list of options.
The outputs are a .log file, containing some information on the
system and the adaptions made to run it in Q, a protein.pdb and
a water.pdb file. The latter two are needed in the next stage.  

-- NOTE --  
The use of molprobity output is currently under construction!  

# 3.setupFEP  
In the next stage we will prepare the input files for the actual
simulations. This stage uses $qligfep/QligFEP.py, you can use
the -h flag to get more specifics on all the input variables.
This folder includes a useful top level script, generate.py, that
takes a .txt file as input that specifies the ligand pairs that
will be used in this simulation. You can simply run:  

    python setup.py

And this should result in a 1.protein and 2.water folder, containing
the two systems necessary to perform the calculation. Simply go into
the folder and run  

    ./FEP_submit.sh

This will either submit your jobs to a slurm queue (see the setup
section in $setupFEP/READE.md for a more detailed description).
Alternatively the jobs can be run on your local machine.  

A typical procedure is to put the two systems in a top layer folder
to easily track running calculations (e.g. 1.TESTRUN) in the case
of this example. Of course, these toplayer scripts can be easily
adjusted for any other use (e.g. to calculate solvation free
energies, you add a water and vacuum leg, instead of protein water.  

# 4.analyzeFEP  
The last step includes running the analyze_allFEP.py script, also
present in the 3.setupFEP folder. In the particular case of this
example, the analysis will be conduced on the 1.TESTRUN folder.
The output is a results.txt file containing the DDG values obtained
for every ligand pair. Exponential averaging (Zwanzig), Overlap
Sampling (OS) and Bennet acceptance ratio (BAR) are used, and the
error is a standard error of the mean over the (finished) replicates.
Alternatively, one can get the free energies per leg by using the
analyze_FEP.py in the main folder. (e.g. you can run:
python $/qligfep/analyze_FEP.py -h to get an overview of the
required input).  
