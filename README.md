# CIPHER: Calphad Integrated PHase-field solvER

## What is this?

CIPHER is a phase-field simulation code for microstructure evolution in multi-component alloy systems. Features include:

- Automatic parallel adaptive mesh refinement.
- Local truccation error estimates and adaptive time stepping.
- Designed for large number of phases (10-10000), with constant computational complexity, i.e. computational cost is independent of the number of phases.
- Efficient grand-canonical-based phase-field implementation with direct use of Compound-Energy-Formalism and other CALPHAD thermodynamic models for multi-component systems [1].
- Designed for MPI parallelization and scalability.

## Requirements

This software requires MPI, p4est v2.2 [2], and PETSc v3.11 [3].

Download the latest p4est source from http://www.p4est.org and unpack, for eg., to ~/software/p4est-2.2. 

Set P4EST_DIR environment variable:
```bash
# replace ~/software/p4est-2.2 with your unpack directory
export P4EST_DIR=~/software/p4est-2.2
```

To install:
```bash
# nagivate to the project directory
cd P4EST_DIR
# configure (on CSF3)
./configure --prefix=$PWD/my-build --enable-mpi CC=mpicc FC=mpif90 F77=mpif77 CXX=mpic++ CFLAGS=-O2 -msse4.2 -axSSE4.2,AVX,CORE-AVX2 CXXFLAGS=-O2 -msse4.2 -axSSE4.2,AVX,CORE-AVX2 FFLAGS=-O2 -msse4.2 -axSSE4.2,AVX,CORE-AVX2 BLAS_LIBS=$MKLROOT/intel64_lin/libmkl_blas95_ilp64.a
# compile
make
make install
```

Download the latest PETSc source from https://www.mcs.anl.gov/petsc/ and unpack, for eg., to ~/software/petsc-3.11.3. 

Set PETSC_DIR and PETSC_ARCH environment variables:
```bash
# replace ~/software/petsc-3.11.3 with your unpack directory, and my-arch with whatever you want to call your installation
export PETSC_DIR=~/software/petsc-3.11.3
export PETSC_ARCH=my-arch
```

To install:

```bash
# nagivate to the project directory
cd PETSC_DIR
# configure (on CSF3)
./configure  --with-pthread --download-metis --download-parmetis --download-chaco --with-mkl_pardiso-dir=$MKLROOT/intel64 --with-mkl_sparse-dir=$MKLROOT/lib/intel64 --with-mkl_sparse_optimize-dir=$MKLROOT/lib/intel64 --download-hypre --download-ml --download-triangle --download-ctetgen --download-hdf5 --with-zlib --with-p4est-dir=$P$EST_DIR/my-build --with-blaslapack-dir=$MKLROOT/lib --with-cxx-dialect=C++11 --with-debugging=0 COPTFLAGS="-O2 -msse4.2 -axSSE4.2,AVX,CORE-AVX2" CXXOPTFLAGS="-O2 -msse4.2 -axSSE4.2,AVX,CORE-AVX2" FOPTFLAGS="-O2 -msse4.2 -axSSE4.2,AVX,CORE-AVX2" PETSC_ARCH=$PETSC_ARCH PETSC_DIR=$PETSC_DIR
# compile
make PETSC_ARCH=$PETSC_ARCH PETSC_DIR=$PETSC_DIR
make install PETSC_ARCH=$PETSC_ARCH PETSC_DIR=$PETSC_DIR
```

## Usage

Make sure PETSC_DIR and PETSC_ARCH have been set correctly (see above). To compile:
```bash
# nagivate to the project directory
cd CIPHER_DIR
make clean
# compile 
make install
```

Running an example:
```bash
# if not already set, add CIPHER_DIR/bin to your PATH
export PATH=CIPHER_DIR/bin:$PATH
# navigate to examples folder
cd examples/GrainBoundaryPrecipitate
# run example
mpiexec -n 4 cipher.exe -geomfile GrainBoundaryPrecipitate.geom -interfacefile N3.interface
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
 
