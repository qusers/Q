# QligFEP solvation free energies test set

This test set contains files for the FEP SER -> ALA perturbation, part of the benchmark set.

Note that for the moment we are testing with 15A spheres

To run this test:
- Untar oldQrunfiles.tar.gz, this contains reference Q files in the old format.
  These files are used as input for the new QligFEP application for backward compatability.

- Change to folder 2.Q-gpu and run setup.py, this will generate all necessary files to
  run the test. Note that you need to add cluster specific properties in the settings.py
  folder

- Submit the jobs with python submit.py
