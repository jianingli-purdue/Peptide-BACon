## PRO-BACon
Scripts to efficiently traverse the peptide search space and find new self-assembling sequences. PRO-BACon includes scripts to build pdb models of peptides and peptide arrays, as well as scripts to simulate these arrays in Desmond and test for aggregation.

PRO-BACon includes a Montecarlo algorithm that searches for self-aggregating peptides of user-defined sequence length via all-atom or all-atom coarse-grain (AACG) simulation.

### Prerequisites
- Schrodinger2018-4
- NumPy

### File Descriptions
- grid.py -- class that represents actual peptide grid/array structure
- box.py -- higher level class that represents a grid structure and its' Desmond simulations
- dev.py -- classes and functions to build individual peptides
- cgSim.py -- class that represents an AACG job
- mdSim.py -- class that represents an all-atom job
- sysBuilder.py -- class that represents the system builder job
- sysMin.py -- script used in Montecarlo algorithm to convert all-atom model to AACG
- montecarlo.py -- Montecarlo algorithm, one sequence at a time
- montecarlo_mix.py -- Montecarlo algorithm, co-assembly of two peptides
- mcplot.py -- functions used to plot SASA/Montecarlo results (useful but you will need to change path variables)
- make_array.py -- script to create pdb of a peptide array
- make_pep.py -- script to create pdb of a single peptide

## Getting Started
### Creating peptide and peptide array pdbs
The make_array.py script accepts either an input file or user-edited variables. Input files must be placed in the 'input/' directory and take two forms. The first uses spacing, height, length, and width. Spacing is the distance between peptides in Angstroms. Height, length, and width are the dimensions of the array in terms of number of peptides - i.e. this array is 30 peptides high. Peptide sequences are given either by sequence and shape (alpha/sheet) or via a pdb structure that is placed in the 'input/' directory.
```
# !Name dimensions_test
# !source aa_ref/
* spacing=15 height=30 length=300 width=300
>
AA alpha
A.pdb
<
```
The second uses values for spacing, total, and height. The total number of peptides divided by the height should have a square root.
```
# !Name dimensions_test
# !source aa_ref/
# !conc 0.8 0.1 0.1
# !rand True
* spacing=15 total=100 height=4
>
AA alpha
AA beta
A.pdb
<
```
Adding blocked=True or rand=True to the * line will return either a blocked or randomized array (note the two are mutually exclusive). If you would like to control concentrations, simply add the concentration to the concentration line. All of these options can also be set by not using an input file and simply editing the box constructor.

NOTE: Concentrations must add to 1, the number of structures you use must be a multiple of the total population, and dimensions must generally make sense!

To make a single peptide, edit the variables in make_pep.py. Both make_array.py and make_pep.py save models to the 'out/' directory. Both files are run by using:
$SCHRODINGER/run make_x.py (inputfile if necessary)

### Running Montecarlo algorithm
Edit the Montecarlo file of your choice to have a path variable pointing to where you would like the jobs stored (make sure this directory exists) and run:
$SCHRODINGER/run montecarlo.py

### Authors
* Marlo Zorman
* Jonathon Ferrell
* Jianing Li
