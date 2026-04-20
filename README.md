# CIPHER: Calphad Integrated PHase-field solvER

## What is this?

CIPHER is a phase-field simulation code for microstructure evolution in multi-component alloy systems. Features include:

- Automatic parallel adaptive mesh refinement.
- Local truncation error estimates and adaptive time stepping.
- Designed for large number of phases (10-10000), with constant computational complexity, i.e. computational cost is independent of the number of phases.
- Efficient grand-canonical-based phase-field implementation with direct use of Compound-Energy-Formalism and other CALPHAD thermodynamic models for multi-component systems [1].
- Designed for MPI parallelization and scalability.

## Requirements

This software requires MPI, p4est [2], and PETSc [3]. 

Note that the following installation instructions are specific to users of the University of Manchester's Computational Shared Facility (CSF).  Non CSF users will need to follow their local procedures to install PETSc.  

To install CIPHER, first clone the repo using

```
git clone --recurse-submodules https://github.com/micmog/CIPHER.git
```

Then `cd` to the `CIPHER` directory you have just cloned, and run the install script:

```bash
cd CIPHER
source ./install.sh
```
Note that you should change the value of `CIPHER_DIR` in the `install.sh` script to the directory you have
cloned this repo into.

A slurm jobscript is also provided to automate the installation at [install-jobscript.sh](install-jobscript.sh).

## Usage

Running an example:
```bash
# Set environment variables
export CIPHER_DIR=<path/to/cipher>
export PATH=$CIPHER_DIR/bin:$PATH
# Navigate to examples folder
cd $CIPHER_DIR/examples/GrainBoundaryPrecipitate
# run example
mpiexec -n 4 cipher.exe --config GrainBoundaryPrecipitate.yaml
```

## Contact

This code is maintained by the Microstructure Modelling Group at the University of Manchester. 
For questions, comments, bug-reports or contributions please email Dr. Pratheek Shanthraj at pratheek.shanthraj@manchester.ac.uk.

## Funding

Development of CIPHER is primarily funded through EPSRC programme grants NEWAM (EP/R027218/1) and LightForm (EP/R001715/1).

## References

[1] Grand-canonical phase-field implementation: https://doi.org/10.1016/j.cma.2020.113029  
[2] p4est: http://www.p4est.org    
[3] PETSc: https://www.mcs.anl.gov/petsc/  

## License

You are free to use and redistribute this is a free software under the terms of the GNU General Public License v3.0.
 
