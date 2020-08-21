This test includes a single point energy calculation for a water sphere of size 20A. 

Reference energy files, as well as initial coordinates and velocities, are taken from a test run with Q5.7, all scripts are in a tarball and can be extracted with:

```bash
tar -xvf testfiles.tar.gz
```

Run the test with the following command:

```bash
qdyn.py -t Q5_data/benzene-water.top -m Q5_data/eq1.inp -d TEST -r Q5_data
```

A new folder called **TEST** should be created after a succesful run, notice the -d flag to qdyn.py, which stands for the name of the folder where results should be created.

To check if the obtained results match the reference results use:

```bash
python3 compare.py
```
