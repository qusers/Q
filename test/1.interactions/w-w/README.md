This test includes a single point energy calculation for a water sphere of size 20A. 

Reference energy files, as well as initial coordinates and velocities, are taken
from a test run with Q5.7, all scripts are in a tarball and can be extracted:

tar -xvf testfiles.tar.gz

Run the test with the following command:

python ../../../bin/qdyn.py -t topo.json -m eq1.inp -d TEST -r -r Q5_data/ --verbose

Then run compare.py to check if the obtained results match the reference results.
