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
To install PETSc:
```bash
#Navigate back to your software folder, then download PETSc with:
cd $HOME/software
git clone https://gitlab.com/petsc/petsc.git petsc
#Now set an environment variable to tell the system where PETSc is located:
export PETSC_DIR=$HOME/software/petsc
export PETSC_ARCH=cipher
#Note that if you put PETSc in a different folder to that used in this example, you will need to alter this variable.
#Load cmake (required to install some external packages):
module load tools/gcc/cmake/3.11.4
#Load MPI compilers:
module load mpi/intel-17.0/openmpi/4.0.1
#Navigate to the folder containing PETSc and configure it, including the necessary flags:
cd $PETSC_DIR
./configure --download-metis --download-parmetis --download-chaco --download-triangle --download-ctetgen --download-pragmatic --download-eigen --download-hypre --download-ml --download-hdf5 --download-zlib --download-yaml --download-p4est --with-pthread --with-mkl_pardiso-dir=$MKLROOT --with-mkl_sparse-dir=$MKLROOT --with-mkl_sparse_optimize-dir=$MKLROOT --with-blaslapack-dir=$MKLROOT --with-cxx-dialect=C++11 --with-debugging=0 COPTFLAGS="-O2 -msse4.2 -axSSE4.2,AVX,CORE-AVX2" CXXOPTFLAGS="-O2 -msse4.2 -axSSE4.2,AVX,CORE-AVX2" FOPTFLAGS="-O2 -msse4.2 -axSSE4.2,AVX,CORE-AVX2" PETSC_ARCH=$PETSC_ARCH PETSC_DIR=$PETSC_DIR
#Then make PETSc:
make PETSC_DIR=$PETSC_DIR PETSC_ARCH=$PETSC_ARCH
make install PETSC_DIR=$PETSC_DIR PETSC_ARCH=$PETSC_ARCH all
```

To install CIPHER:
```bash
#Now navigate to your software folder and download CIPHER with:
cd $HOME/software
git https://github.com/micmog/CIPHER.git
#Now put CIPHER in your path, so that it can be called from anywhere in your system.
export CIPHER_DIR=$HOME/software/CIPHER
PATH=$CIPHER_DIR/bin:$PATH:$HOME/bin
# nagivate to the project directory
cd $CIPHER_DIR
# compile 
make clean
make install
```

Setting up environment to use CIPHER:
```bash
export OMP_NUM_THREADS=1
module load tools/env/proxy
module load mpi/intel-17.0/openmpi/4.0.1
export PETSC_DIR=$HOME/software/petsc
export PETSC_ARCH=cipher
export CIPHER_DIR=$HOME/software/CIPHER
PATH=$PETSC_DIR/$PETSC_ARCH/bin:$CIPHER_DIR/bin:$PATH
LD_LIBRARY_PATH=$PETSC_DIR/$PETSC_ARCH/lib:$LD_LIBRARY_PATH
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
 
