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

To install CIPHER:
```bash
#Set location to install CIPHER
export CIPHER_DIR=<path/to/cipher>
#Replace <path/to/cipher> with specific location, e.g. export CIPHER_DIR=$HOME/software/CIPHER. It is recommended Place this in your .bashrc/.bash_profile.
#Clone the CIPHER repository with:
git clone --recurse-submodules https://github.com/micmog/CIPHER.git $CIPHER_DIR
#Install PETSc:
cd $CIPHER_DIR/petsc
module load tools/gcc/cmake/3.11.4
module load mpi/intel-17.0/openmpi/4.0.1
./configure --download-metis --download-parmetis --download-chaco --download-triangle --download-ctetgen --download-pragmatic --download-eigen --download-hypre --download-ml --download-hdf5 --download-zlib --download-yaml --download-p4est --with-pthread --with-mkl_pardiso-dir=$MKLROOT --with-mkl_sparse-dir=$MKLROOT --with-mkl_sparse_optimize-dir=$MKLROOT --with-blaslapack-dir=$MKLROOT --with-cxx-dialect=C++11 --with-debugging=0 COPTFLAGS="-O2 -msse4.2 -axSSE4.2,AVX,CORE-AVX2" CXXOPTFLAGS="-O2 -msse4.2 -axSSE4.2,AVX,CORE-AVX2" FOPTFLAGS="-O2 -msse4.2 -axSSE4.2,AVX,CORE-AVX2" PETSC_ARCH=cipher PETSC_DIR=$CIPHER_DIR/petsc
make PETSC_ARCH=cipher PETSC_DIR=$CIPHER_DIR/petsc
make install PETSC_ARCH=cipher PETSC_DIR=$CIPHER_DIR/petsc all
cd $CIPHER_DIR
#Set up CIPHER environment
source $CIPHER_DIR/load_CIPHER.sh
#Install CIPHER 
make clean
make install
```

## Usage

Running an example:
```bash
# after setting up environment to use CIPHER, navigate to examples folder
cd $CIPHER_DIR/examples/GrainBoundaryPrecipitate
# run example
mpiexec -n 4 cipher.exe --config GrainBoundaryPrecipitate.yaml
```

## Contact

This code is maintained by the Microstructure Modelling Group at the University of Manchester. 
For questions, comments, bug-reports or contributions please email Dr. Pratheek Shanthraj at pratheek.shanthraj@manchester.ac.uk.

## References

[1] Grand-canonical phase-field implementation: https://arxiv.org/abs/1906.10503  
[2] p4est: http://www.p4est.org    
[3] PETSc: https://www.mcs.anl.gov/petsc/  

## License

You are free to use and redistribute this is a free software under the terms of the GNU General Public License v3.0.
 
