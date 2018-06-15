# Qfep version history

## Release 5.7, 2016-01-07
- Johan Sund's code for computing free energies using Thermodynamic
  Integration, Overlap Sampling, and Bennet's Acceptance Ratio (BAR)
  has been included.

## Release 4.00, 2000-02-08
- Is identical to 3.75

## Release 3.75, 2000-02-08
- Average energies printed now reflect the skipping of data points. 
- Previously, the averages were always over all the data in the energy
  file.
- Improved the handling of read errors due to wrong number of
  off-diagonal functions in the energy file.

## Release 3.4, 1999-04-15
- Dynamic allocation of all large arrays avoids excessive memory usage
  and  makes the program run faster (less swapping).

## Release 3.32, 1999-04-14
- Added 'Bin-averaged summary'.

## Release 3.31, 1999-04-13
- Removed hard-coded maximum number of energy files (=100). 

## Release 3.3, 1999-03-31
- Changed format of output to be able to use the output as is in
  gnuplot:  

    gnuplot index 0 = average energies (one long line for each state in each file)
	gnuplot index 1 = FEP summary
	gnuplot index 2 = reaction free energy

Writes prompt to standard error if possible. This means one can run:  

    Qfep3 > fep.log

and still run the program interactively.

## Release 3.2, 1998-05-09
- Changed format of off-diagonals in energy file to be same as old Qdyn/Qfep.
- Now reports module versions.
