# CIPHER: Calphad Integrated PHase-field solvER

## What is this?

CIPHER is a phase-field simulation code for microstructure evolution in multi-component alloy systems. Features include:

- Automatic parallel adaptive mesh refinement.
- Local truncation error estimates and adaptive time stepping.
- Designed for large number of phases (10-10000), with constant computational complexity, i.e. computational cost is independent of the number of phases.
- Efficient grand-canonical-based phase-field implementation with direct use of Compound-Energy-Formalism and other CALPHAD thermodynamic models for multi-component systems [1].
- Designed for MPI parallelization and scalability.

## Requirements

This software requires MPI, p4est v2.2 [2], and PETSc v3.12 [3]. 

Note that the following installation instructions are specific to users of the University of Manchester's Computational Shared Facility (CSF).  Non CSF users will need to follow their local procedures to install, p4est and PETSc.  
To install p4est:
```bash
#Load MPI compilers and OpenBLAS with the command:
module load mpi/intel-17.0/openmpi/3.1.3
module load libs/gcc/openblas/0.3.6
#Navigate to your software folder.  If you don’t have one enter
mkdir $HOME/software
cd $HOME/software
#Download p4est with the command line utility wget:
wget http://p4est.github.io/release/p4est-2.2.tar.gz
#The p4est download will be compressed in a tar file.  Extract this using:  
tar -xvf p4est-2.2.tar.gz
#You can now delete the p4est tar file:
rm p4est-2.2.tar.gz
#Those without the wget utility can download the latest p4est source from http://www.p4est.org and unpack, for eg., to ~/software/p4est-2.2.
#Navigate into the extracted p4est folder.  Now configure the p4est software library, with the correct command arguments (flags):
./configure --prefix=$PWD/intel-17.0-openblas --enable-mpi CC=mpicc FC=mpif90 F77=mpif77 CXX=mpic++ CFLAGS='-O2 -msse4.2 -axSSE4.2,AVX,CORE-AVX2' CXXFLAGS='-O2 -msse4.2 -axSSE4.2,AVX,CORE-AVX2' FFLAGS='-O2 -msse4.2 -axSSE4.2,AVX,CORE-AVX2' BLAS_LIBS=$OPENBLASDIR/lib/libopenblas.a
#This runs a script that will localize the source distribution so that it will compile and load on your local system.  Next enter:
make
make install
#Now you need to add an environment variable to your system where p4est is located
export P4EST_DIR=$HOME/software/p4est-2.2/intel-17.0-openblas 
```

To install PETSc:
```bash
#Navigate back to your software folder, then download PETSc with:
cd $HOME/software
wget http://ftp.mcs.anl.gov/pub/petsc/release-snapshots/petsc-lite-3.12.1.tar.gz
#and extract it with the tar utility:
tar -xvf petsc-lite-3.12.1.tar.gz
#Now set an environment variable to tell the system where PETSc is located:
export PETSC_DIR=$HOME/software/petsc-3.12.1
export PETSC_ARCH=intel-17.0-mkl
#Note that if you put PETSc in a different folder to that used in this example, you will need to alter this variable.
#Make a folder, inside your PETSc folder and name if after the architecture that you will be using (eg mkdir intel-17.0-mkl).  Make another folder, inside the architecture folder for your external packages and name it ‘external packages’ (eg mkdir intel-17.0-mkl/externalpackages).  Add another line to .bash_profile to set an environment variable to specify what architecture PETSc will use:
mkdir $PETSC_DIR/$PETSC_ARCH/externalpackages
cd $PETSC_DIR/$PETSC_ARCH/externalpackages
#You need to download three external packages for PETSc, so navigate into the external packages folder.  The first of these is Triangle:
wget http://ftp.mcs.anl.gov/pub/petsc/externalpackages/Triangle.tar.gz
#next one is HDF5
wget https://support.hdfgroup.org/ftp/HDF5/prev-releases/hdf5-1.8/hdf5-1.8.18/src/hdf5-1.8.18.tar.gz
#next one is YAML
wget http://pyyaml.org/download/libyaml/yaml-0.1.4.tar.gz
#and finally Chaco:
wget  http://ftp.mcs.anl.gov/pub/petsc/externalpackages/Chaco-2.2-p2.tar.gz
#Load cmake (required to install some external packages):
module load tools/gcc/cmake/3.11.4
#Navigate to the folder containing PETSc and configure it, including the necessary flags:
cd $PETSC_DIR
./configure  --with-pthread --download-yaml=$PWD/$PETSC_ARCH/externalpackages/yaml-0.1.4.tar.gz --download-metis --download-parmetis --download-chaco=$PWD/$PETSC_ARCH/externalpackages/Chaco-2.2-p2.tar.gz --with-mkl_pardiso-dir=$MKLROOT --with-mkl_sparse-dir=$MKLROOT --with-mkl_sparse_optimize-dir=$MKLROOT --download-hypre --download-ml --download-triangle=$PWD/$PETSC_ARCH/externalpackages/Triangle.tar.gz --download-ctetgen --download-hdf5=$PWD/$PETSC_ARCH/externalpackages/hdf5-1.8.18.tar.gz --with-zlib --with-p4est-dir=$P4EST_DIR --with-blaslapack-dir=$MKLROOT --with-cxx-dialect=C++11 --with-debugging=0 COPTFLAGS="-O2 -msse4.2 -axSSE4.2,AVX,CORE-AVX2" CXXOPTFLAGS="-O2 -msse4.2 -axSSE4.2,AVX,CORE-AVX2" FOPTFLAGS="-O2 -msse4.2 -axSSE4.2,AVX,CORE-AVX2" PETSC_ARCH=$PETSC_ARCH PETSC_DIR=$PETSC_DIR
#Then make PETSc:
make PETSC_DIR=$PETSC_DIR PETSC_ARCH=$PETSC_ARCH
make install PETSC_DIR=$PETSC_DIR PETSC_ARCH=$PETSC_ARCH all
```

To install CIPHER:
```bash
#To allow git to function, load the tools/env/proxy module:
module load tools/env/proxy
#Now navigate to your software folder and make a folder for CIPHER (mkdir CIPHER).  Go to the CIPHER github page (https://github.com/micmog/CIPHER).  Click on the button ‘Clone or Download’ and select ‘Download ZIP’.  Go to your Downloads folder, move the download to your software folder and rename it CIPHER.  You can also gain CIPHER by navigating to your software folder and using Git clone:
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

Setting up environment:
```bash
# .bash_profile is a file that runs every time you open a bash shell. It is convenient to define all the environment variables there so it's available everytime you log on to the CSF.
# Open .bash_profile for editing, with the following command:
vim $HOME/.bash_profile
#Then add this line to the file:
module load tools/env/proxy
module load mpi/intel-17.0/openmpi/3.1.3
export P4EST_DIR=$HOME/software/p4est-2.2/intel-17.0-openblas 
export PETSC_DIR=$HOME/software/petsc-3.12.1
export PETSC_ARCH=intel-17.0-mkl
export CIPHER_DIR=$HOME/software/CIPHER
PATH=$CIPHER_DIR/bin:$PATH
```

## Usage

Running an example:
```bash
# if not already set, add CIPHER_DIR/bin to your PATH
export PATH=CIPHER_DIR/bin:$PATH
# navigate to examples folder
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
 
