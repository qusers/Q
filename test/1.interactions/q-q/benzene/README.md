This test includes a single point energy calculation for a water sphere of size 20A. 

Reference energy files, as well as initial coordinates and velocities, are taken
from a test run with Q5.7, all scripts are in a tarball and can be extracted:

tar -xvf testfiles.tar.gz

Run the test with the following command:

python ../../../../bin/qdyn.py -t benzene-vacuum.json -m eq1.json -d TEST -r Q5_data -f FEP1.json/

Then to check if the obtained results match the reference results:

python compare.py
